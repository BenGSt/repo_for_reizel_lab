#!/usr/bin/env bash


main() #<sample_dir>
{
  sample_dir=$1
  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
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

  cd "$sample_dir" || exit 1
  time deduplicate_bismark ./*bismark*bam
  rm $(find . -name '*.bam' | grep -v deduplicated) #delete bam with duplicates

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

main "$@"