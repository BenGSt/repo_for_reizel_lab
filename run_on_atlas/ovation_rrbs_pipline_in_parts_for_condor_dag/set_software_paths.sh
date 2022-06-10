#!/bin/bash

set_software_paths()
{
  #	echo setting software paths
	# on Atlas I'm using conda for env managment, so the software should be in path already.
	# leaving this function if I want to change something without the env.
	#TODO: consider using calls to program name from the rest of the script and delete this function.

	PYTHON3=python

	TRIM_GALORE=trim_galore

	BOWTIE2=bowtie2

	SAMTOOLS=samtools

	BISMARK=bismark

	FASTQC=fastqc

	DIVERSITY_TRIM_SCRIPT=/Local/bfe_reizel/anaconda3/envs/ovation_rrbs_pipeline_2022/bin/trimRRBSdiversityAdaptCustomers.py

	BEDTOOLS=bedtools

	#tecan nudup tool for pcr duplicates
	NUDUP=/Local/bfe_reizel/anaconda3/envs/ovation_rrbs_pipeline_2022/bin/nudup.py
}

