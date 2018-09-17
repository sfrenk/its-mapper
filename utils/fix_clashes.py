#!/usr/bin/env python3

### Resolve clashes (ITS sites called on both strands in overlapping locations) ###

import argparse
import sys

parser = argparse.ArgumentParser(description = "Resolve clashes (ITS sites called on both strands in overlapping locations)")

parser.add_argument("file", help = "bed file containing ITS sites")
parser.add_argument("-o", "--output", help = "output filename", default = "score_merged_final.bed")
parser.add_argument("-s", "--score_col", help = "column containing score (default: 5", default = 5, type = int)

args = parser.parse_args()

outfile = open(args.output, "w")
logfile = open(args.output + ".log", "w")
score_index = args.score_col - 1

with open(args.file) as f:

	prev_line = f.readline()

	for line in f:

		# Get previous coords
		prev_chrom = prev_line.strip().split("\t")[0]
		prev_start = int(prev_line.strip().split("\t")[1])
		prev_end = int(prev_line.strip().split("\t")[2])
		prev_strand = prev_line.strip().split("\t")[3]
		prev_score = float(prev_line.strip().split("\t")[score_index])

		chrom = line.strip().split("\t")[0]
		start = int(line.strip().split("\t")[1])
		end = int(line.strip().split("\t")[2])
		strand = line.strip().split("\t")[3]
		score = float(line.strip().split("\t")[score_index])
		
		# Check for clash with previous line
		if (chrom != prev_chrom or start >= prev_end):
			# No clash
			outfile.write(prev_line)
			prev_line = line
		
		else:
			# Clash
			logfile.write(prev_line + line)
			
			if score > prev_score:
				# Current ITS wins.
				prev_line = line

			elif score == prev_score:
				# Tie
				print("ERROR: exactly equal scores! Not sure what to do here...")
				print(prev_line)
				print(line)
				sys.exit(1)


outfile.close()
logfile.close()
