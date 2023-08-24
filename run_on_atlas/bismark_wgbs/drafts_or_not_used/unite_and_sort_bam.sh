#!/usr/bin/env bash

# This script is used to unite and sort bam files from bismark alignment
# Author: Ben Steinberg
# Date: 22/08/2023

cd $samp_dir || exit 1
sample_name=$(basename $samp_dir)
merged_bam=${sample_name}_bismark_bt2.bam
#merge all the bam files
samtools merge -f $(find ./split -name "*.bam") $merged_bam
samtools sort -n $merged_bam -o ${merged_bam/.bam/_name_sorted.bam}
