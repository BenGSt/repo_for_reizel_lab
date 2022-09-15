#!/bin/bash
input_file=$1
output_file=$2
GTF=/storage/bfe_reizel/bengst/genomic_reference_data/hg38/hg38.p13/genes/hg38.knownGene.gtf

source /Local/bfe_reizel/anaconda3/bin/activate rna-seq_hisat2_htseq_deseq2_2022
htseq-count --stranded no $input_file $GTF > $output_file