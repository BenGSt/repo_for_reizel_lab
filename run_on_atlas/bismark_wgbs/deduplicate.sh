#!/usr/bin/env bash


main() #<sample_dir> <split>
{
  sample_dir=$1
  split=$2 #USAGE: set if input are multiple split fastq files (if run_data_split.sh was run)
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
  if [[ $split ]]; then
    deduplicate_bismark --multiple $(ls ./*bismark*bam | sort)
  else
    deduplicate_bismark ./*bismark*bam
  fi
  rm -v $(find . -name '*.bam' | grep -v deduplicated) #delete bam file(s) with duplicates

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