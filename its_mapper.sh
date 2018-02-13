#!/usr/bin/bash -e
#SBATCH -n 4
#SBATCH --array=0-102
#SBATCH -o logs/%a.out
#SBATCH -e logs/%a.err
#SBATCH -t 1-0

module add python bedtools samtools emboss r

#### Set up directories ####

directories=("scoring_window" "hmm" "bg")

for directory in ${directories[@]}; do
	if [[ ! -d $directory ]]; then
		mkdir $directory
	fi
done

#line_number=$((${SLURM_ARRAY_TASK_ID} + 1))
#chromosome=$(sed -n "${line_number}p" ../chunks.bed | cut -f1)

#### Scoring window ####

scoring_window -w 100 -c chunks.bed -o scoring_window/score_${SLURM_ARRAY_TASK_ID}_ -i $SLURM_ARRAY_TASK_ID -x

#### HMM ####

Rscript ~/scripts/its_mapper/call_its_sites.R -o hmm/hmm_${SLURM_ARRAY_TASK_ID}_fwd.txt scoring_window/score_${SLURM_ARRAY_TASK_ID}_fwd.txt

Rscript ~/scripts/its_mapper/call_its_sites.R -o hmm/hmm_${SLURM_ARRAY_TASK_ID}_rev.txt scoring_window/score_${SLURM_ARRAY_TASK_ID}_rev.txt

#### Define ITS blocks ####
strands=("fwd" "rev")

for i in ${strands[@]}; do
	# Merge consequetive coords into regions
	bedtools merge -i hmm/hmm_${SLURM_ARRAY_TASK_ID}_${i}.txt > hmm/hmm_${SLURM_ARRAY_TASK_ID}_${i}.bed
	# Assign score to each region
	bedtools map \
	-a hmm/hmm_${SLURM_ARRAY_TASK_ID}_${i}.bed \
	-b scoring_window/score_${SLURM_ARRAY_TASK_ID}_${i}.txt \
	-c 4 \
	-o mean > bg/${SLURM_ARRAY_TASK_ID}_${i}.bg
done
	