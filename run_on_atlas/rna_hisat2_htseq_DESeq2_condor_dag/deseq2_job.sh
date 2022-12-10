#!/bin/bash

# copy all htseq-count ouput files to one dir
mkdir ./htseq_output
cp `find -name *htseq-count_output.txt` ./htseq_output

source /Local/bfe_reizel/anaconda3/bin/activate rna_hisat2_htseq_deseq2_2022_v2
Rscript $REPO_FOR_REIZEL_LAB/scripts_2022/deg_deseq2.R "$@"
