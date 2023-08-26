#!/bin/bash

N_CORES=10
N_PARALLEL_INSTANCES=4 #each instance uses ~3 cores
BUFFER_SIZE=4G #buffer size for unix sort
MEM=4GB


help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh, align_to_genome.sh
	resources: $N_CORES cores, $MEM RAM

  produces coverage report (only CpG for to include all Cytosines use option --CX)

	-output-dir
	#TODO: -keep-files #OT OB *.bedGraph.gz


	optional:
  -extra-options "multiple double quoted options" (e.g. --ignore <int>. See Bismark manual)
  -bam-dir <path> Use for m-bias fix (when the bam files are in a different directory than the output directory).
EOF
}


main()
{
  arg_parse "$@"
	source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
	cd "$output_dir" || exit 1
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

#	time methylation_calling
	#TODO try without time , see if it speeds up or less memory used
  methylation_calling

	#cleanup
  rm -v $(find ./ | grep -P 'OT|OB')
  rm -v $(find . -name "*.bedGraph.gz")

	echo
	echo
	echo \################################
	echo \################################
	echo finished: $script_name
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
  #if bam_dir is an empty string (i.e. not set), find the bam files in the current directory
  if [[ -z $bam_dir ]]; then
    alignment_output=$(find . -name '*bismark*deduplicated*bam')
  else
    alignment_output=$(find $bam_dir -name '*bismark*deduplicated*bam')
  fi
  echo $alignment_output | grep 'pe' && paired="-p" || paired=""
  command=$(echo bismark_methylation_extractor --bedgraph $paired $ignore_r2 --multicore $N_PARALLEL_INSTANCES --gzip --buffer_size $BUFFER_SIZE $extra $alignment_output)

  #--ample_memory speeds things up for samples over 10 million reads or so. since it may take over an hour to get going ATLAS policy holds the jobs.
#  command=$(echo bismark_methylation_extractor --ample_memory --bedgraph $paired $ignore_r2 --multicore $N_PARALLEL_INSTANCES --gzip  $extra $alignment_output)
  echo $SCRIPT_NAME runnig: $command
	$command
}


arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -h|--help)
        help
        exit 1
        ;;
      -output-dir)
        output_dir="$2"
        shift
        shift
        ;;
      -ignore_r2)
        ignore_r2=$(echo --ignore_r2 "$2")
        shift
        shift
        ;;
      -extra-options)
        extra=$2
        shift
        shift
        ;;
      -bam-dir)
        bam_dir=$2
        shift
        shift
        ;;
      *)
        echo ERROR: unknown option: $1
        help
        exit 1
        ;;
    esac
  done
}


main "$@"