#!/usr/bin/env python3

import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(description = "Make bed file containing chunk coordinates for sliding window")

parser.add_argument("-s", "--chunk_size", help = "Size of genome chunk (default: 1000000 bases)", type = int, default = 1000000)
parser.add_argument("-o", "--output", help = "Output filename (Default: chunks.bed", default = "chunks.bed")

args = parser.parse_args()

# Get chromosome size info

#chrom_sizes_file = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.bedtools_genome_file"

chrom_sizes_file = "/home/sfrenk/Documents/Resources/Seq/WS251/bedtools_genome_file.txt"

chrom_sizes = OrderedDict()

with open(chrom_sizes_file) as f:
	for line in f:
		chrom_sizes[line.strip().split()[0]] = int(line.strip().split()[1])

# Write chunks to bed file
outfile = open(args.output, "w")

for chrom, max_coord in chrom_sizes.items():
	
	# 0-based index
	start_coord = 0
	end_coord = args.chunk_size
	while end_coord < max_coord:
		outfile.write("\t".join([chrom, str(start_coord), str(end_coord)]) + "\n")
		start_coord = start_coord + args.chunk_size
		end_coord = end_coord + args.chunk_size

	outfile.write("\t".join([chrom, str(start_coord), str(max_coord)]) + "\n")

outfile.close()