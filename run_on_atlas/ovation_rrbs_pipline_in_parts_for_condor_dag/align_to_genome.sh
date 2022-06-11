#!/bin/bash

GENOMIC_REFERENCE_LOCATION=/storage/bfe_reizel/bengst/genomic_reference_data
BISMARK_GENOME_LOCATION=${GENOMIC_REFERENCE_LOCATION}/from_huji/mm10/Sequence/WholeGenomeFasta


help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh
		resources: 20 cores, 40GB RAM

	-single-end or -paired-end
	-input_fastq_file <sample.fq.gz> or -paired_input_fastq_files <sample_R1.fq.gz> <sample_R2.fq.gz>
	-n_cores
EOF
}


main()
{
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

	arg_parse "$@"
	set_software_paths

	time align_to_genome

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


align_to_genome()
{
  #see http://felixkrueger.github.io/Bismark/Docs/ :
    #"--parallel 4 for e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~48GB of RAM,
    # but at the same time reduce the alignment time to ~25-30%. You have been warned."
  n_parallel_instances = $(( $n_cores / 5 ))

  if [[ $read_type == "single_end" ]] ; then
    trim_diversity_output=$(echo $trim_galore_output | sed 's/\.gz/_trimmed.fq.gz/')
    command=$(echo $BISMARK --multicore $n_parallel_instances --bowtie2 $BISMARK_GENOME_LOCATION $trim_diversity_output)
	else
    trim_diversity_output_1=$(echo $trim_galore_output_1 | sed 's/\.gz/_trimmed.fq.gz/')
    trim_diversity_output_2=$(echo $trim_galore_output_2 | sed 's/\.gz/_trimmed.fq.gz/')
    command=$(echo $BISMARK --multicore $n_parallel_instances --bowtie2 $BISMARK_GENOME_LOCATION -1 $trim_diversity_output_1 -2 $trim_diversity_output_2)
	fi

  echo runnig: $command \($(date)\)
  $command

	#ASK_TZACHI: Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)
	#is this what we want?
}


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
	-n_cores)
        n_cores="$2"
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