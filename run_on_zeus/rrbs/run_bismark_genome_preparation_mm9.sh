#!/bin/bash

#PBS  -N  bismark_genome_preperation
#PBS  -q  zeus_all_q
#PBS  -m  abe
#PBS  -M  s.benjamin@technion.ac.il
#PBS  -l select=1:ncpus=20
#PBS  -l select=mem=64gb
#PBS  -l walltime=24:00:00

PBS_O_WORKDIR=/home/s.benjamin/genomic_reference_data/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/
cd $PBS_O_WORKDIR

eval "$(micromamba shell hook --shell=bash)"
micromamba activate  /home/s.benjamin/micromamba/envs/wgbs_bismark_pipeline_2023
echo running"bismark_genome_preparation --bowtie2 --parallel 10 --verbose /home/s.benjamin/genomic_reference_data/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ >> bismark_genome_preparation.log 2>&1" > bismark_genome_preparation.log
bismark_genome_preparation --bowtie2 --parallel 10 --verbose /home/s.benjamin/genomic_reference_data/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ >> bismark_genome_preparation.log 2>&1