#!/usr/bin/env python3

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
import subprocess

parser = argparse.ArgumentParser(description = "Get telomere homology scores for ITS sequences")

parser.add_argument("file", help = "Bed file containing ITS sites")
parser.add_argument("-o", "--output_file", help = "Bed file containing ITS sites annotated with homology score", default = "its_homology.bed")
parser.add_argument("-r", "--ref", help = "Reference fasta", required = True)
parser.add_argument("-t", "--telo_repeat", help = "Telomere repeat sequence (default: TTAGGG)", default = "TTAGGG")
parser.add_argument("-s", "--size_filter", help = "Minimum size for ITS (default: 0)", type = int, default = 0)


args = parser.parse_args()

# Get reverse complement of telo sequence
telo_seq = args.telo_repeat
comp = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
telo_comp = [ comp[x] for x in telo_seq ]
telo_revcomp = "".join(telo_comp)[::-1]


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
			seq = subprocess.check_output(["samtools", "faidx", args.ref, chrom+ ":" + str(start) + "-" + str(end)]).decode("utf-8")

			# Get telomere ref sequence
			ref_length = int(math.ceil(float(size/float(6))))
			if strand == "+":
				telo_ref = telo_seq * ref_length
			elif strand == "-":
				telo_ref = telo_revcomp * ref_length
			else:
				raise Exception("ERROR: strand must be + or -")

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
