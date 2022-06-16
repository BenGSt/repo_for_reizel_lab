#!/bin/bash

N_CORES=10
N_PARALLEL_INSTANCES=4 #each instance uses 3 cores (according to doc) #ignores n_cores
BUFFER_SIZE=10G #buffer size for unix sort
MEM=10GB


help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh, align_to_genome.sh
	resources: $N_N_CORES cores, $MEM RAM

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

	time methylation_calling

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


methylation_calling()
{
  if [[ $read_type == "single_end" ]]; then
    trim_galore_output=$(echo $input_fastq |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_trimmed.fq.gz/')
    trim_diversity_output=$(echo $trim_galore_output | sed 's/\.gz/_trimmed.fq.gz/')
    rename=$(echo $trim_diversity_output| sed 's/\.fq_trimmed/_trimmed/')
	  alignment_output=$(echo $rename | sed 's/\.fq\.gz/_bismark_bt2.bam/')
	  #By default, this mode will only consider cytosines in CpG context, but it can be extended to cytosines in any sequence context by using the option --CX
    command=$(echo bismark_methylation_extractor --multicore $N_PARALLEL_INSTANCES --bedGraph --buffer_size $BUFFER_SIZE --output methylation_extractor_output $alignment_output)
	else
	  trim_galore_output_1=$(echo $input_fastq_1 |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_val_1.fq.gz/')
    rename=$(echo $trim_diversity_output| sed 's/\.fq_trimmed/_trimmed/')
    trim_diversity_output_1=$(echo $trim_galore_output_1 | sed 's/\.gz/_trimmed.fq.gz/')
    rename_1=$(echo $trim_diversity_output_1| sed 's/\.fq_trimmed/_trimmed/')
	  alignment_output=$(echo $rename_1 | sed 's/\.fq\.gz/_bismark_bt2_pe.bam/')
#	  mv $alignment_output $(echo $alignment_output | sed 's/_R1_001_val_1//')
#	  alignment_output=$(echo $alignment_output | sed 's/_R1_001_val_1//')
	  # TODO: looks like bismark wrote one bam for the 2 paired files,
	    # but it's name includes the R1 name I changed the name - make sure the bam file really represents both
    command=$(echo bismark_methylation_extractor -p --multicore $N_PARALLEL_INSTANCES --bedGraph --buffer_size $BUFFER_SIZE --output methylation_extractor_output $alignment_output)
	fi

	echo $SCRIPT_NAME runnig: $command
	$command
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