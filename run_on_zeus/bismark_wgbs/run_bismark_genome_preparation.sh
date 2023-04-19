#!/bin/bash

#PBS  -N  bismark_genome_preperation
#PBS  -q  zeus_all_q
#PBS  -m  abe
#PBS  -M  s.benjamin@technion.ac.il
#PBS  -l select=1:ncpus=20
#PBS  -l select=mem=64gb
#PBS  -l walltime=24:00:00
PBS_O_WORKDIR=/home/s.benjamin/extra_storage_1TB/genomic_reference_data/hg38
cd $PBS_O_WORKDIR

conda activate /home/s.benjamin/micromamba/envs/wgbs_bismark_pipeline_2023
bismark_genome_preparation --bowtie2 --parallel 10 --verbose /home/s.benjamin/extra_storage_1TB/genomic_reference_data/hg38