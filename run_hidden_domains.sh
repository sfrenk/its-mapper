#!/usr/bin/bash
#SBATCH --mem 30G
#SBATCH -t 1-0

module add r perl

# Convert bedGraph to hiddenEnrichment input file
# + strand
awk -v OFS="\t" '{print NR,$1,$3,$4}' /nas/longleaf/home/sfrenk/pine/scoring_window/its_window_scores.bedGraph | sed '1iID\tchr\tpos\tcount' > its_window_bins_fwd.txt

hiddenEnrichment -g /nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.bedtools_genome_file -b its_window_bins_fwd.txt -o its_results_fwd

# - strand
awk -v OFS="\t" '{print NR,$1,$3,$5}' /nas/longleaf/home/sfrenk/pine/scoring_window/its_window_scores.bedGraph | sed '1iID\tchr\tpos\tcount' > its_window_bins_rev.txt

hiddenEnrichment -g /nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.bedtools_genome_file -b its_window_bins_rev.txt -o its_results_rev
