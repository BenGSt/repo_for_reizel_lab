#!/bin/bash

N_CORES=1
MEM=500MB

help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh
	resources: 1 core, 500MB RAM

	-single-end or -paired-end
	-input_fastq_file <sample.fq.gz> or -paired_input_fastq_files <sample_R1.fq.gz> <sample_R2.fq.gz>
	-output_dir
EOF
}


main()
{
  arg_parse "$@"
	set_software_paths
	cd $output_dir
	script_name=$(echo $0 | awk -F / '{print $NF}')

	echo
	echo
	echo \################################
	echo \################################
	echo running: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo

	time trim_diversity_adaptors

	echo
	echo
	echo \################################
	echo \################################
	echo finished: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo
}


trim_diversity_adaptors()
{
	#https://github.com/nugentechnologies/NuMetRRBS#diversity-trimming-and-filtering-with-nugens-diversity-trimming-scripts
	if [[ $read_type == "single_end" ]] ; then
	  trim_galore_output=$(echo $input_fastq |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_trimmed.fq.gz/')
	  echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output
	  python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output
	else
	  trim_galore_output_1=$(echo $input_fastq_1 |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_val_1.fq.gz/')
	  trim_galore_output_2=$(echo $input_fastq_2 |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_val_2.fq.gz/')
	  echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output_1 -2 $trim_galore_output_2
	  python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output_1 -2 $trim_galore_output_2
	fi
}


set_software_paths()
{
  #	echo setting software paths
	# on Atlas I'm using conda for env managment, so the software should be in path already.
	# leaving this function if I want to change something without the env.
	#TODO: consider using calls to program name from the rest of the script and delete this function.

  source /Local/bfe_reizel/anaconda3/bin/activate ovation_rrbs_pipeline_2022

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



arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
     -single-end)
        read_type="single_end"
        shift # past argument
        ;;
     -paired-end)
        read_type="paired_end"
        shift # past argument
        ;;
     -input_fastq_file)
        input_fastq="$2"
        shift # past argument
        shift # past value
        ;;
      -paired_input_fastq_files)
        input_fastq_1="$2"
        shift # past argument
        input_fastq_2="$2"
        shift # past argument2
        shift # past value
        ;;
	-output_dir)
        output_dir="$2"
        shift # past argument
        shift # past value
        ;;
      -*|--*)
        help
        exit 1
        ;;
      -h|--help)
        help
        exit 1
        ;;
    esac
  done
}


main "$@"