#!/bin/bash

# copy all htseq-count ouput files to one dir
cp `find -name *htseq-count_output.txt` ./htseq_output

source /Local/bfe_reizel/anaconda3/bin/activate rna_hisat2_htseq_deseq2_2022_v2
Rscript /srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/run_deseq2.R "$@"
