#!/bin/bash
input_file=$1
output_file=$2
genome=$3

if [[ $genome == "mm10" ]]; then
  gtf=/storage/bfe_reizel/bengst/genomic_reference_data/mm10.ensGene.gtf
elif [[ $genome == "hg38" ]]; then
  gtf=/storage/bfe_reizel/bengst/genomic_reference_data/hg38/hg38.p13/genes/hg38.ensGene.gtf
else
  echo genome not recognized
  exit 1
fi

echo running $0 "$@"
echo input_file=$1
echo output_file=$2
echo genome=$3
echo gtf=$gtf


source /Local/bfe_reizel/anaconda3/bin/activate rna-seq_hisat2_htseq_deseq2_2022
htseq-count --stranded no $input_file $gtf > $output_file