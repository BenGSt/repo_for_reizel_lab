#!/usr/bin/env bash

# This script writes and submits submission files for preparing the sample prep jobs, one job per sample is run in order
# to count the reads in the raw data (fastq.gz) file, and write the appropriate submission files needed in order to split
# the input and run the pipeline.

# Once the prep jobs have finished running, another script is used to write the submission file for multiqc, and the top
# level dag files, and submit the WGBS bismark pipline jobs.

# After the WGBS pipeline jobs have finished running, another script is used to fix m-bias by repeating the pipeline,
# starting from the methylation calling stage, and ignoring a specified number of base pairs from read ends.

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab
source $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/shared.sh
source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023

prepapre_sample() {
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

write_prep_submission_files() {
  for sample_name in $(find -L $raw_data_dir -type d | awk -F / 'NR>1{print $NF}' | sort); do
    cat <<EOF >condor_submission_files/prep/${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/prepare_jobs.sh
Arguments = "$@" -job -sample-name $sample_name
request_cpus = 1
RequestMemory = 500MB
universe = vanilla
log = $(pwd)/logs/multiqc_job.log
output = $(pwd)/logs/multiqc_job.out
error = $(pwd)/logs/multiqc_job.out
queue

EOF
  done
}

main() {
  mkdir -p logs
  mkdir -p condor_submission_files/prep/
  n_reads_per_chunk=100000000 #default value (may be overwritten by arg_parse)
  arg_parse "$@"
  if [[ $job -eq 1 ]]; then
    prepare_sample
  else
    echo DEBUG: save cmd
    save_cmd prep.cmd
    write_prep_submission_files "$@"
  fi
}

arg_parse() {
  if [[ $# -eq 0 ]]; then
    help
    exit 1
  fi
  while [[ $# -gt 0 ]]; do
    case $1 in
    -h | --help)
      help
      exit 1
      ;;
    -single-end)
      single_end=1
      shift
      ;;
    -paired-end)
      single_end=0
      shift
      ;;
    -non-directional)
      non_directional="--non_directional"
      shift
      ;;
    -raw-data-dir)
      raw_data_dir=$2
      shift
      shift
      ;;
    -keep-bam)
      keep_bam="-keep-bam"
      shift
      ;;
    -keep-trimmed-fq)
      keep_trimmed_fq="-keep-trimmed-fq"
      shift
      ;;
    -dovetail) #seems this is on by default.
      dovetail="-dovetail"
      shift
      ;;
    -genome)
      genome=$2
      shift
      shift
      ;;
    -n-reads-per-chunk) #for splitting fastq files, default is 100M
      n_reads_per_chunk=$2
      shift
      shift
      ;;
    -extra-trim-galore-options)
      extra_trim_opts=$(echo -extra-trim-galore-options \'"$2"\')
      shift
      shift
      ;;
    -extra-meth_extract-options)
      extra_meth_opts=$(echo -extra-options \'"$2"\')
      shift
      shift
      ;;
    -job)
      job=1
      shift
      ;;
    -sample-name)
      sample_name=$2
      shift
      shift
      ;;
    *)
      help
      exit 1
      ;;
    esac
  done
}

main "$@"
