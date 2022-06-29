#!/bin/bash

N_CORES=10
MEM=16GB
N_PARALLEL_INSTANCES=2


GENOMIC_REFERENCE_LOCATION=/storage/bfe_reizel/bengst/genomic_reference_data
BISMARK_GENOME_LOCATION=${GENOMIC_REFERENCE_LOCATION}/from_huji/mm10/Sequence/WholeGenomeFasta


help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh
		resources: 10 cores, 16GB RAM

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

	time align_to_genome
	cleanup

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
  # Atlas max cpu request is 10 so I want to have 2 instances of bismark (5 cores each theoretically)
  # This is set in align_jobs.sub .


  if [[ $read_type == "single_end" ]] ; then
    trim_galore_output=$(echo $input_fastq |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_trimmed.fq.gz/')
    trim_diversity_output=$(echo $trim_galore_output | sed 's/\.gz/_trimmed.fq.gz/')
    rename=$(echo $trim_diversity_output| sed 's/\.fq_trimmed/_trimmed/')
    mv $trim_diversity_output $rename

    command=$(echo $BISMARK --multicore $N_PARALLEL_INSTANCES --bowtie2 $BISMARK_GENOME_LOCATION $rename)
	else
	  trim_galore_output_1=$(echo $input_fastq_1 |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_val_1.fq.gz/')
	  trim_galore_output_2=$(echo $input_fastq_2 |awk -F / '{print $NF}'| sed 's/\(\.fastq\|.fq\)\.gz/_val_2.fq.gz/')
    trim_diversity_output_1=$(echo $trim_galore_output_1 | sed 's/\.gz/_trimmed.fq.gz/')
    trim_diversity_output_2=$(echo $trim_galore_output_2 | sed 's/\.gz/_trimmed.fq.gz/')
    rename_1=$(echo $trim_diversity_output_1| sed 's/\.fq_trimmed/_trimmed/')
    rename_2=$(echo $trim_diversity_output_2| sed 's/\.fq_trimmed/_trimmed/')
    mv $trim_diversity_output_1 $rename_1
    mv $trim_diversity_output_2 $rename_2

    command=$(echo $BISMARK --multicore $N_PARALLEL_INSTANCES --bowtie2 $BISMARK_GENOME_LOCATION -1 $rename_1 -2 $rename_2)
	fi

  echo runnig: $command \($(date)\)
  $command

	#ASK_TZACHI: Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)
	#is this what we want?
}


cleanup()
{
  cmd="rm $rename $rename_1 $rename_2"
  echo cleanup: "$cmd"
  $cmd

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