#!/usr/bin/env bash

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab
source $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/shared.sh --source-only
source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023

prepapre_sample() {
  raw_dir=$1
  sample_name=$2

  mkdir -p condor_submission_files/$sample_name
  mkdir -p logs/$sample_name

  # count reads to see if fastq file is longer than n_reads_per_chunk reads, if so it will be split into chunks of
  # length n_reads_per_chunk, given as an argument to this script, defaults to 100M.
  count_reads

  #write sub files for trimming and aligning each chunk (or one sub file if no splitting)
  write_trim_and_align_sub_files

  # the following sub files are not dependent on splitting
  write_deduplicate_job_submission_file
  write_methylation_calling_job_submission_file
  write_bam2nuc_job_submission_file
  write_make_tiles_job_submission_file
  write_sample_dag_file
}

# This script is writes the assortment of submission files needed to run one sample as a job on Atlas.
# Another script is used to write the submission file for multiqc,and top level dag files.
# Yet another script provides the user cli for preparing jobs for all the samples.
main() {
  n_reads_per_chunk=100000000 #default value (may be overwritten by arg_parse)
  prepapre_sample "$@"
}


main "$@"
