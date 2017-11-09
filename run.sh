#!/usr/bin/bash
#SBATCH --mem 4G
#SBATCH -t 0-5
#SBATCH --array 0-102
#SBATCH -o logs/%a.out
#SBATCH -e logs/%a.err

module add python samtools emboss
scoring_window -c chunks.bed -i $SLURM_ARRAY_TASK_ID -x

cat score_* | sort -k 1,1 -k 2,2n > its_window_scores.bedGraph
