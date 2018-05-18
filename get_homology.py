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

parser = argparse.ArgumentParser(description = "Get telomere homology scores for ITS sequences")

parser.add_argument("file", help = "Bed file containing ITS sites")
parser.add_argument("-o", "--output_file", help = "Bed file containing ITS sites annotated with homology score", default = "its_homology.bed")
parser.add_argument("-s", "--size_filter", help = "Minimum size for ITS (default: 0)", type = int, default = 0)
parser.add_argument("-sp", "--species", help = "Species (elegans (default), briggsae or remani)", default = "elegans")


args = parser.parse_args()
outfile = open(args.output_file, "w")

with open(args.file) as f:
	for line in f:
		chrom = line.strip().split("\t")[0]
		start = int(line.strip().split("\t")[1]) + 1
		end = int(line.strip().split("\t")[2])
		strand = line.strip().split("\t")[3]
		size = end-start

		if size >= args.size_filter:
			# Get ITS sequence
			seq = samtools_lookup.get_seq(chrom, start, end, name = None, zero = False, species = args.species).seq

			# Get telomere ref sequence
			ref_length = int(math.ceil(float(size/float(6))))
			if strand == "+":
				telo_ref = "TTAGGC" * ref_length
			elif strand == "-":
				telo_ref = "GCCTAA" * ref_length
			else:
				print("ERROR: strand must be + or -")
				sys.exit(1)

			# Perform alignment with water
			with open("its_seq.temp", "w") as fi:
				fi.write(str(seq))
			with open("telo.temp", "w") as ft:
				ft.write(telo_ref)

			water_cmd = WaterCommandline(gapopen = 10, gapextend = 0.5, asequence = "its_seq.temp", bsequence = "telo.temp", stdout = True, auto = True)
			stdout, stderr = water_cmd()
			identity = re.findall("# Identity:.*\((.+)\%\)", stdout)[0]

			outfile.write(line.strip() + "\t" + str(identity) +"\n")

outfile.close()
os.remove("its_seq.temp")
os.remove("telo.temp")
