library(depmixS4)
library(tidyverse)
library(argparse)

parser <- ArgumentParser(description = "Get a bed file of its sites from scoring_window.py output")

parser$add_argument("file", help = "bedgraph file output of scores for a single chromosome")
parser$add_argument("-o", "--output", help = "output filename (default: hmm.bed)", default = "hmm.bed")

args <- parser$parse_args()

# Read in data
data <- read.table(args$file)
colnames(data) <- c("chrom", "start", "end", "signal")

# Only one chromosome please
chromosome <- unique(data$chrom)
if (length(chromosome) != 1){
    print("ERROR: Input file should contain data for one chromosome")
    q(save = "no", status = 1)
}

# Build model
mod <- depmix(signal ~ 1, data=data, nstates=2)

# Fit model
f <- fit(mod)
summary(f)

# Get states
esttrans <- posterior(f)

# make bed file containing all positions with state = 2
data_out_bed <- cbind(data, esttrans) %>% filter(state == 2) %>% dplyr::select(chrom, start, end)
write.table(data_out_bed, args$output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    