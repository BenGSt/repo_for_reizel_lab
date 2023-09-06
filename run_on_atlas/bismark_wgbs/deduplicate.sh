#!/usr/bin/env bash
#  ----------------------------------------
#  Project: Reizel Lab Bioinformatics Pipelines
#  Pipeline: Bismark WGBS
#  Script: deduplicate.sh
#  Author: Ben G. Steinberg
#  Last Update: 4 Sep 2023
#  ----------------------------------------

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

# USAGE: deduplicate.sh <sample_dir> {--paired-end or --single-end} <split>
#   sample_dir: directory containing the alignment output bam file(s)
#   --paired-end or --single-end: specify if the input data is paired end or single end
#   split: set if input are multiple split fastq files (if split_fastq.sh was run)
main() { #<sample_dir> {--paired-end or --single-end} <split>
  sample_dir=$1
  split=$3 #USAGE: set if input are multiple split fastq files (if run_data_split.sh was run)
  script_name=$(echo $0 | awk -F / '{print $NF}')
  print_info "running: " "$script_name " "$@"

  if [[ $2 == "-paired-end" ]]; then
    flags="-p"
  elif [[ $2 == "-single-end" ]]; then
    flags="-s"
  else
    echo "ERROR: must specify -paired-end or -single-end"
    exit 1
  fi

  if [[ $split ]]; then
    flags=$(echo $flags "--multiple")
  fi

  cd "$sample_dir" || exit 1

  echo "If output files from previous runs exist, they will be removed as to not corrupt the current run."
  rm -fv $(find . -name "*deduplicated*bam")

  # NOTE: with the wgbs_bismark_pipeline_2023 conda environment, we get broken pipe errors from samtools and perl.
  # This is a known issue and should not affect the output. The samtools error can be fixed by using an older version
  # of samtools (tested with  -samtools_path /Local/bfe_reizel/samtools-0.1.19/). The perl error might also be fixed by
  # using an older version of perl, but this has not been tested.
  deduplicate_bismark $flags $(find . -name "*bismark*bam" | sort) || exit 1

  # rename deduplicated bam file to remove chunk and val from name of paired end
  # and chunk and trimmed from se (for multiqc)
  dedup_bam=$(ls *deduplicated*bam)
  new_bam_name=${dedup_bam/_chunk_[0-9]*_val_[0-9]/}
  new_bam_name=${new_bam_name/_chunk_[0-9]*_trimmed/}
  new_bam_name=${new_bam_name/_R[1-2]/}
  mv $dedup_bam $new_bam_name

  rm -v $(find . -name '*.bam' | grep -v deduplicated) #delete bam file(s) with duplicates

  print_info "finished: " "$script_name " "$@"
}

main "$@"
