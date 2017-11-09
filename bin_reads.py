#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description = "Bin data in bedGraph file")

parser.add_argument("file", help = "Input bedGraph file")
parser.add_argument("-b", "--binsize", help = "Size of bins in bases (default: 1000", type = int, default = 1000)
parser.add_argument("-o", "--output", help = "Output filename (default: binned.bg)", default = "binned.bg")

args = parser.parse_args()

outfile = open(args.output, "w")

with open(args.file) as f:

	prev_chrom = ""
	bin_start = 0
	prev_start = 0
	prev_end = 0
	bin_score = 0

	for line in f:
		
		chrom = line.strip().split("\t")[0]
		start = int(line.strip().split("\t")[1])
		end = int(line.strip().split("\t")[2])
		score = float(line.strip().split("\t")[3])

		if chrom != prev_chrom or end-bin_start > args.binsize:

			if prev_chrom != "":
				
				# End of bin
				# Output current bin
				mean_count = bin_score/(prev_end-bin_start)
				outfile.write("\t".join([chrom, str(bin_start), str(prev_end), str(mean_count)]) + "\n")

				# New bin
				if chrom != prev_chrom:

					# Reached the end of the chromosome
					bin_start = 0
				
				else:

					bin_start = prev_end
				
				bin_score = 0

		prev_chrom = chrom
		prev_start = start
		prev_end = end
		bin_score = bin_score + score

	# Output any remaining piece of genome at the end of the last chromosome
	mean_count = bin_score/(prev_end-bin_start)
	outfile.write("\t".join([chrom, str(bin_start), str(prev_end), str(mean_count)]) + "\n")

outfile.close()
