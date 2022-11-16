#!/bin/bash
r1=$1
r2=$2
output_file=$3
summary_file=$4 # {}_hisat2.summary.txt
#human
hisat2_idx=/storage/bfe_reizel/bengst/genomic_reference_data/hg38/analysisSet/chromosomes_hisat2_index/hisat2-index

#mouse
#hisat2_idx=/storage/bfe_reizel/bengst/genomic_reference_data/from_huji/mm10/Sequence/whole_genome_fasta_hisat2_index/hisat2-index

source /Local/bfe_reizel/anaconda3/bin/activate rna-seq_hisat2_htseq_deseq2_2022
mkdir -p $(dirname $summary_file)
hisat2 -p 10 -x $hisat2_idx -1 $r1 -2 $r2 --summary-file $summary_file | samtools sort -n --output-fmt BAM > $output_file
#TODO: single end read option