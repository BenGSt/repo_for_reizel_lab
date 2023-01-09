#!/bin/bash

N_CORES=10
N_PARALLEL_INSTANCES=4 #each instance uses ~3 cores
BUFFER_SIZE=10G #buffer size for unix sort
MEM=3GB


help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh, align_to_genome.sh
	resources: $N_N_CORES cores, $MEM RAM

	-output-dir
	-keep-bam
	-keep-trimmed-fq
EOF
}


main()
{
  arg_parse "$@"
	source /Local/bfe_reizel/anaconda3/bin/activate ovation_rrbs_pipeline_2022
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

	time methylation_calling
  cleanup

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
  alignment_output=$(find . -name '*bismark*deduplicated*bam')
  echo $alignment_output #debug
  echo $alignment_output | grep 'pe' && paired="-p --ignore_r2 2" || paired=""
  command=$(echo bismark_methylation_extractor $paired --multicore $N_PARALLEL_INSTANCES --gzip --bedGraph --buffer_size $BUFFER_SIZE --output methylation_extractor_output $alignment_output)
  echo $SCRIPT_NAME runnig: $command
	$command
}


cleanup()
{
  rm_bam="rm $alignment_output"
  rm_OT_OB="rm $(find ./ | grep -P 'OT|OB')"
  rm_fq="rm *.fq" #the non gz trimmed fq
  if [[ $keep_bam -eq 0 ]]; then
    $rm_bam
  fi

  if [[ $keep_trimmed_fq -eq 0 ]]; then
    $rm_fq
  fi

  $rm_OT_OB
}


arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -output-dir)
        output_dir="$2"
        shift
        shift
        ;;
      -keep-bam)
        keep_bam=1
        shift
        ;;
      -keep-trimmed-fq)
        keep_trimmed_fq=1
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