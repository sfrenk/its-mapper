#!/usr/bin/env python3

import argparse
from collections import OrderedDict
import sys

parser = argparse.ArgumentParser(description = "Make bed file containing chunk coordinates for sliding window")

parser.add_argument("-s", "--chunk_size", help = "Size of genome chunk (default: 1000000 bases)", type = int, default = 1000000)
parser.add_argument("-w", "--window_size", help = "window for sliding_window.py. The chunk size must be divisible by the window size (default = 100", type = int, default = 100)
parser.add_argument("-g", "--genome_file", help = "bedtools genome file", default = "/home/sfrenk/Documents/Resources/Seq/WS251/bedtools_genome_file.txt")
parser.add_argument("-o", "--output", help = "Output filename (Default: chunks.bed", default = "chunks.bed")

args = parser.parse_args()

# Check chunk/window size
if args.chunk_size % args.window_size != 0:
	print("ERROR: Chunk size must be divisible by window size")
	sys.exit(1)

# Get chromosome size info
chrom_sizes = OrderedDict()

with open(args.genome_file) as f:
	for line in f:
		chrom_sizes[line.strip().split()[0]] = int(line.strip().split()[1])

# Write chunks to bed file
outfile = open(args.output, "w")

for chrom, max_coord in chrom_sizes.items():
	
	# 0-based index
	start_coord = 0
	end_coord = start_coord + args.chunk_size + args.window_size - 1
	while end_coord < max_coord:
		outfile.write("\t".join([chrom, str(start_coord), str(end_coord)]) + "\n")
		start_coord = start_coord + args.chunk_size
		end_coord = end_coord + args.chunk_size

	# When the distance between the chunk stat and the end of the chromosome is less than the chunk size, output the remaining region
	outfile.write("\t".join([chrom, str(start_coord), str(max_coord)]) + "\n")

outfile.close()