#!/usr/bin/env python3

import argparse
import re
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("file", help = "genome fasta file")
parser.add_argument("-o", "--output", help = "output filename", default = None)

args = parser.parse_args()

# Setup input file handle
# Input file may or may not be gzipped
if args.file.endswith(".gz"):
	f = gzip.open(args.file, "rb")
else:
	f = open(args.file)

# Set up output file handle
if args.output is None:
	output_filename = str(re.search("([^/]+)\.fa(.gz)?$", args.file).group(1)) + ".chromosome_sizes.txt"
else:
	output_filename = args.output
outfile = open(output_filename, "w")

name = None
count = 0

for line in f:

	try:
		# Line needs to be decoded if input is gzipped
		line = line.decode("utf-8")
	except:
		pass
	
	if line.startswith(">"):
		if name is not None:
			outfile.write(name + "\t" + str(count) + "\n")
		name = line.strip()[1:]
		count = 0
	else:
		count = count + len(line.strip())

f.close()

# Output final sequence
if count != 0:
	outfile.write(name + "\t" + str(count) + "\n")

outfile.close()
