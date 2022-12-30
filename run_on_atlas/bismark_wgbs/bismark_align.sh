#!/bin/bash

N_CORES=10
MEM=16GB
N_PARALLEL_INSTANCES=2 # ~5 cores per instance



help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh
		resources: $N_CORES cores, $MEM RAM

	<-single-end> or <-paired-end>
	<-output-dir>
	-genome <mm10 or hg38>
	[-non-directional]   instructs Bismark to use all four alignment outputs (OT, CTOT, OB, CTOB)
EOF
}


main()
{
  arg_parse "$@"
	cd "$output_dir" || exit 1
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

	if [[ $3 == "mm10" ]]; then
	  bismark_genome_location=/storage/bfe_reizel/bengst/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta
	elif [[ $3 == "hg38" ]]; then
	  bismark_genome_location=/storage/bfe_reizel/bengst/genomic_reference_data/hg38/analysisSet/hg38.analysisSet.chroms/
	else
	  echo genome not recognized
	  exit 1
  fi

	time align_to_genome

	echo
	echo
	echo \#################################
	echo \#################################
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
    trim_galore_output=$(find . -name '*trimmed.fq*')
    command=$(echo bismark --multicore $N_PARALLEL_INSTANCES --bowtie2 $bismark_genome_location $trim_galore_output $non_directional)
	else
	  trim_galore_output_1=$(find . -name '*val_1.fq*')
	  trim_galore_output_2=$(find . -name '*val_2.fq*')
    command=$(echo bismark --multicore $N_PARALLEL_INSTANCES --bowtie2 $bismark_genome_location -1 $trim_galore_output_1 -2 $trim_galore_output_2 $non_directional)
	fi

  echo runnig: $command
  $command
}


cleanup()
{
  cmd="rm $rename $rename_1 $rename_2"
  echo cleanup: "$cmd"
  $cmd

}

arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -single-end)
        read_type="single_end"
        shift
        ;;
      -paired-end)
        read_type="paired_end"
        shift
        ;;
      -output-dir)
        output_dir="$2"
        shift
        shift
        ;;
      -non-directional)
        non_directional="--non_directional"
        shift
        ;;
      -genome)
        genome=$2
        shift
        shift
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