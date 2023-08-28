#!/usr/bin/env bash


main() #<sample_dir> <split> {--paired-end|--single-end}
{
  sample_dir=$1
  split=$2 #USAGE: set if input are multiple split fastq files (if run_data_split.sh was run)
#TODO: meth extractor is broken, commenting this out for now to see if it works without it
#  if [[ $3 == "-paired-end" ]]; then
#    pe_or_se="-p"
#  elif [[ $3 == "-single-end" ]]; then
#    pe_or_se="-s"
#  else
#    echo "ERROR: must specify -paired-end or -single-end"
#    exit 1
#  fi

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
    deduplicate_bismark --multiple $(find . -name "*bismark*bam" | sort) $pe_or_se
  else
    deduplicate_bismark ./*bismark*bam $pe_or_se
  fi

  # rename deduplicated bam file to remove chunk and val from name of paired end
  # and chunk and trimmed from se (for multiqc)
  dedup_bam=$(ls *deduplicated*bam)
  new_bam_name=${dedup_bam/_chunk_[0-9]*_val_[0-9]/}
  new_bam_name=${new_bam_name/_chunk_[0-9]*_trimmed/}
  new_bam_name=${new_bam_name/_R[1-2]/}
  mv $dedup_bam $new_bam_name

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
