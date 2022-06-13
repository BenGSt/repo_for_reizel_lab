#!/bin/bash

N_CORES=8
MEM=300MB

help()
{
	cat << EOF
	run first
	resources: 8 cores, 300MB RAM

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

	if [[ $read_type == "single_end" ]]
	then
		time trim_illumina_adapter_single_end $input_fastq
	else #if [[ $read_type == "paired_end" ]]
		time trim_illumina_adapter_paired_end $input_fastq_1 $input_fastq_2
	fi

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


trim_illumina_adapter_paired_end() #<R1> <R2>
{
	#positional argument are  R1, R2 fastq file to trim
	#note on multicores from cutadapt manual
		#To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
	#note from trim_galore manual
		#It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
		#--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
	echo runnig: trim_galore --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $1 $2 --cores $N_CORES --fastqc \($(date)\)
	${TRIM_GALORE} --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $1 $2 --cores $N_CORES --fastqc
}


trim_illumina_adapter_single_end()
{
	#first positional argument is fastq file to trim
	#note on nulticores from cutadapt manual
		#To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
	#note from trim_galore manual
		#It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
		#--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
	echo runnig: trim_galore --adapter AGATCGGAAGAGC $1 --cores $N_CORES \($(date)\)
	${TRIM_GALORE} --adapter AGATCGGAAGAGC $1  --cores $N_CORES --fastqc
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