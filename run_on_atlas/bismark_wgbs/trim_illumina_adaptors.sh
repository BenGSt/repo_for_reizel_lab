#!/bin/bash

N_CORES=8
MEM=300MB

help()
{
	cat << EOF
	run first
	resources: $N_CORES cores, $MEM RAM

	-single-end or -paired-end
	-input_fastq_file <sample.fq.gz> or -paired_input_fastq_files <sample_R1.fq.gz> <sample_R2.fq.gz>
	-output_dir
EOF
}


main()
{
  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
	arg_parse "$@"
	mkdir -p $output_dir
  cd $output_dir
	script_name=$(echo $0 | awk -F / '{print $NF}')

	echo
	echo
	echo \#################################
	echo \#################################
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
	echo \#################################
	echo \#################################
	echo finished: $script_name "$@"
	echo date: $(date)
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
	cmd="trim_galore --dont_gzip --paired $1 $2 --cores $N_CORES --fastqc"
	echo runnig: $cmd
	$cmd
}


trim_illumina_adapter_single_end() #<R1>
{
  cmd="trim_galore $1 --cores $N_CORES --fastqc"
	echo runnig: $cmd
	$cmd
}





arg_parse()
{
  if [[ $# -eq 0 ]]; then
    help
    exit 1
  fi
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
	   -output-dir)
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