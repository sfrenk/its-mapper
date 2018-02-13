#!/usr/bin/env python3

import samtools_lookup
import argparse
import math
import sys
# pairwise2 currently gives segfault error
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
from Bio.Emboss.Applications import WaterCommandline
from Bio import AlignIO
import re
import os

parser = argparse.ArgumentParser(description = "Sliding window score of genome chunk")

parser.add_argument("-c", "--chunk_file", help = "Bed file containing coordinates of genomic chunks", required = True)
parser.add_argument("-i", "--chunk_index", help = "Genomic chunk index", type = int, required = True)
parser.add_argument("-w", "--window_size", help = "Size of sliding window (default: 6 bases", type = int, default = 6)
parser.add_argument("-o", "--output", help = "Output file prefix (Default: score_<index>_ where <index> is the chunk index", default = "DEFAULT")
parser.add_argument("-x", "--clean_mode", help = "Remove all temporary files", default = False, action = "store_true")
parser.add_argument("-g", "--gapopen", help = "Gap open penalty (default: 10)", type = int, default = 10)
parser.add_argument("-e", "--gapextend", help = "Gap extend penalty (default: 0.5)", type = int, default = 0.5)

args = parser.parse_args()


### CHUNK RETREIVAL ###

# Get coordinates of chunk

with open(args.chunk_file) as f:
	chunks = [line for line in f.readlines()]

chunk = chunks[args.chunk_index]

chrom = chunk.strip().split()[0]
start = chunk.strip().split()[1]
end = chunk.strip().split()[2]

# Fetch chunk sequence from the genome fasta file

seq = samtools_lookup.get_seq(chrom, start, end, name = None, zero = True).seq


### OUTPUT FILENAME ###

if args.output == "DEFAULT":
	outfile_prefix = "score_" + str(args.chunk_index) + "_"
else:
	outfile_prefix = args.output


### TELO REF PREPARATION ###

# Make reference sequence - telomere repeat sequence the same length as the query.
# Have to do both telomere strands (TTAGGC/GCCTAA)

telo_repeats = []
for telo_repeat in ["TTAGGC", "GCCTAA"]:

	ref_length = int(math.ceil(float(args.window_size/float(6))))
	ref_seq = telo_repeat * ref_length
	ref_seq = ref_seq[0:args.window_size]

	telo_repeats.append(ref_seq)

with open("telo_fwd_" + str(args.chunk_index) + ".txt", "w") as f:
	f.write(telo_repeats[0] + "\n")

with open("telo_rev_" + str(args.chunk_index) + ".txt", "w") as f:
	f.write(telo_repeats[1] + "\n")


### SLIDING WINDOW PREPARATION ###

# Sliding window - Make multifasta file containing windowed sequence
# Also make lists of genomic regions for output later

with open("query_" + str(args.chunk_index) + ".fa", "w") as f:
	current_pos = 0
	starts = []
	ends = []

	while current_pos + (args.window_size) <= len(seq):
		
		# Calculate the genome position (start of chunk + position with chunk)
		genome_pos = int(start) + current_pos

		# The score is going to be assigned to the middle base of the window
		starts.append(genome_pos + round(args.window_size/2))
		ends.append(genome_pos + 1 + round(args.window_size/2))
		query = seq[current_pos:current_pos + args.window_size]
		f.write(">" + str(current_pos) + "\n" + str(query) + "\n")
		current_pos = current_pos + 1

	# Output the last window of sequence
	#starts.append(int(start) + current_pos + round(args.window_size/2))

	# Assign the same score to each base within the entire window. This shouldn't matter if chunk size is divisible by window size. If the chunk is at the end of a chromosome, the final window will be perfect telomere anyway
	#ends.append(int(end))
	#query = seq[current_pos:]
	#f.write(">" + str(current_pos) + "\n" + str(query) + "\n")


### ALIGNMENT ###

# Alignment chunk sequences to perfect telomere with water algorithm
# Need to align with both forward and reverse telomere strand

telo_scores = {}

for telo in ["telo_fwd_", "telo_rev_"]:

	water_cmd = WaterCommandline(gapopen = args.gapopen, gapextend = args.gapextend, asequence = telo + str(args.chunk_index) + ".txt", bsequence = "query_" + str(args.chunk_index) + ".fa", stdout = True, auto = True)
	stdout, stderr = water_cmd()

	# Extract scores from water output

	##
	with open("water_res" + telo + ".txt", "w") as fy:
		fy.write(stdout)
	##

	scores = re.findall("# Score: ([0-9.]+)", stdout)
	telo_scores[telo] = scores

### OUTPUT ###

# fwd
with open(outfile_prefix + "fwd.txt", "w") as f:
	for i in range(len(scores)):
		f.write("\t".join([chrom, str(starts[i]), str(ends[i]), str(telo_scores["telo_fwd_"][i])]) + "\n")

# rev
with open(outfile_prefix + "rev.txt", "w") as f:
    for i in range(len(scores)):
        f.write("\t".join([chrom, str(starts[i]), str(ends[i]), str(telo_scores["telo_rev_"][i])]) + "\n")

### CLEANUP ###

if args.clean_mode:
	for tempfile in ["telo_fwd_" + str(args.chunk_index) + ".txt", "telo_rev_" + str(args.chunk_index) + ".txt", "query_" + str(args.chunk_index) + ".fa"]:
		os.remove(tempfile)
