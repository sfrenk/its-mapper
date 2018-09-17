#!/usr/bin/bash
set -eu -o pipefail

usage="
    Locate interstitial telomere sites (ITSs) in a genome. The output file, its.bed, contains four columns:

    1. Chromosome
    2. Start (0-indexed)
    3. End
    4. Strand
    5. Number of telomeric repeats
    6. Telomere homology score

    USAGE
       
       bash its_mapper.sh [options]  

    ARGUMENTS
        -g/--genome_fasta
        Uncompressed genome fasta file (required)

        -t/--telo_sequence
        Telomere repeat sequence (default: TTAGGG)

        -b/--bowtie_index
        Bowtie index created from fasta (optional)
    "


# Defaults
script_dir="$(dirname $0)"
utils_dir=${script_dir}/utils
bowtie_index="/proj/seq/data/hg38_UCSC/Sequence/BowtieIndex/genome"
#genome_fasta="/home/sfrenk/Documents/Resources/Seq/human/hg38/genome.fa.gz"
genome_fasta="/proj/ahmedlab/steve/seq/human/hg38/genome.fa"
telo_repeat=TTAGGG

# Parse command line parameters

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -g|--genome_fasta)
        genome_fasta="$2"
        shift
        ;;
        -t|--telo_repeat)
		telo_repeat="$2"
		shift
		;;
		-b|--bowtie_index)
		bowtie_index="$2"
		shift
		;;
    esac
    shift
done


dirs=("bowtie" "ref" "temp")

for i in ${dirs[@]}; do
	if [[ ! -d "$i" ]]; then
		mkdir "$i"
	fi
done


if [[ bowtie_index == "" ]]; then
	# Build bowtie index, genome index and genome file
	echo "Building bowtie index..."
	gunzip -c $genome_fasta > ref/genome.fa
	bowtie-build ref/genome.fa bowtie/genome
	bowtie_index="bowtie/genome"
fi

echo "Indexing genome file"
#samtools faidx $genome_fasta

echo "Building genome file..."
#python3 ${utils_dir}/make_genome_file.py -o ref/genome_file.txt $genome_fasta

# find all telo repeats
echo "Mapping telomere repeats..."
bowtie -S -c -v 0 -a "$bowtie_index" "$telo_repeat" temp/telo_single.sam

samtools view -bh temp/telo_single.sam | samtools sort -o temp/telo_single.bam -
samtools index temp/telo_single.bam

# Merge telo repeats that are whithin 100bp of each other
echo "Creating merged bed file..."
bedtools merge -s -d 100 -i temp/telo_single.bam > temp/single_telo_repeat_merged.bed

# Count the number of telomere repeats in each block
awk -v OFS="\t" '{print $1,$2,$3,".",".",$4}' temp/single_telo_repeat_merged.bed | bedtools intersect -a - -b temp/telo_single.bam -c -s | cut -f 1,2,3,6,7 > temp/single_telo_repeat_merged_count.bed 

# Retain sites with at least 4 hexamers
awk '$5>=4' temp/single_telo_repeat_merged_count.bed > temp/single_telo_repeat_merged_count4.bed

# Remove telomeres and fix clashes
echo "Processing ITS sites..."
awk -v OFS="\t" '{print $1,0,50"\n"$1,$2-50,$2}' ref/genome_file.txt > ref/potential_telomeres.bed
bedtools intersect -v -a temp/single_telo_repeat_merged_count4.bed -b ref/potential_telomeres.bed > temp/its_notelomere.bed
python3 ${utils_dir}/fix_clashes.py -o temp/its_noclash.bed temp/its_notelomere.bed

# Get homology scores
echo "Obtaining telomere homology scores..."
python3 ${utils_dir}/get_homology.py -o its.bed -r $genome_fasta -t $telo_repeat temp/its_notelomere.bed
