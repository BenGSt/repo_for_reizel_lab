#!/usr/bin/env bash
# USAGE: run_fix_mbias.sh --biased_dir <path> --ignore_r1 <int> --ignore_r2  --ignore_3prime <int> --ignore_3prime_r2 <int> [--output-dir <output_dir>]

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab
source $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/run_dag_per_sample_split.sh --source-only

main() {
  bias_fix=1 #for write_sample_dag_file()
  arg_parse "$@"
  extra_meth_opts="$ignore_r1 $ignore_r2 $ignore_3prime $ignore_3prime_r2"
  if [[ ! $output_dir ]]; then
    echo output_dir=${biased_dir}_mbias_fixed
    output_dir=${biased_dir}_mbias_fixed
  fi
  mkdir -p $output_dir
  cd $output_dir
  mkdir -p logs


  #write bismark_methylation_extractor sub files
    #find the sample directories (the first node in the path) for which bam files exist
    for sample_name in $(find $biased_dir -name "*deduplicated*bam" | awk -F / '{print $(NF-1)}'); do
    {
      echo DEBUG: sample_name=$sample_name
      echo
      unset split sep chunk
      sample_names+=($sample_name)
      mkdir -p condor_submission_files/$sample_name
      mkdir -p logs/$sample_name

      # the following sub files are not dependent on splitting
      write_methylation_calling_job_submission_file
      write_bam2nuc_job_submission_file
      write_make_tiles_job_submission_file
      write_sample_dag_file
    }
  done

  write_top_level_dag

  #list jobs and the commands to run them
  #ask if user wants to run them now, if so, run them.
  cat << EOF

Unless you need them, it is recommended to delete the bam files when you are done.

To do so, run: rm -v $(find $biased_dir -name "*deduplicated*bam")

Please download your data and delete it from atlas as soon as you are done.
!Good luck and happy clustering!
EOF
}

help() {
  cat << EOF
  run_fix_mbias.sh --biased_dir <path> {at least one of: --ignore_r1 <int> --ignore_r2  --ignore_3prime <int> --ignore_3prime_r2 <int>} [--output-dir <output_dir>]
obligatory options:
   --biased_dir <path>
at least one of the following:
   --ignore <int>
   --ignore_3prime <int>

   --ignore_r2 <int>
   --ignore_3prime_r2 <int>
non-obligatory options:
   [--output-dir <output_dir>]

options to ignore edges of reads (from Bismark manual):
=====================================
--ignore <int>
    Ignore the first <int> bp from the 5' end of Read 1 (or single-end alignment files) when processing
    the methylation call string. This can remove e.g. a restriction enzyme site at the start of each read or any other
    source of bias (such as PBAT-Seq data).

--ignore_r2 <int>
    Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing results only. Since the first couple of
    bases in Read 2 of BS-Seq experiments show a severe bias towards non-methylation as a result of end-repairing
    sonicated fragments with unmethylated cytosines (see M-bias plot), it is recommended that the first couple of
    bp of Read 2 are removed before starting downstream analysis. Please see the section on M-bias plots in the Bismark
    User Guide for more details.

--ignore_3prime <int>
    Ignore the last <int> bp from the 3' end of Read 1 (or single-end alignment files) when processing the methylation
    call string. This can remove unwanted biases from the end of reads.

--ignore_3prime_r2 <int>
    Ignore the last <int> bp from the 3' end of Read 2 of paired-end sequencing results only. This can remove unwanted
    biases from the end of reads.
EOF
}

arg_parse() {
  if [[ $# -eq 0 ]]; then
    help
    exit 1
  fi
  while [[ $# -gt 0 ]]; do
    case $1 in
      --biased_dir)
        biased_dir="$(realpath $2)"
        shift
        shift
        ;;
      --ignore)
        ignore_r1="--ignore $2"
        shift
        shift
        ;;
      --ignore_r2)
        ignore_r2="--ignore_r2 $2"
        shift
        shift
        ;;
      --ignore_3prime)
        ignore_3prime="--ignore_3prime $2"
        shift
        shift
        ;;
      --ignore_3prime_r2)
        ignore_3prime_r2="--ignore_3prime_r2 $2"
        shift
        shift
        ;;
      --output-dir)
        output_dir=$2
        shift
        shift
        ;;
      *)
        echo "Unknown option: $1"
        help
        exit 1
        ;;
    esac
  done

  if [[ ! $biased_dir ]]; then
    echo "--biased_dir is obligatory"
    help
    exit 1
  fi
}

main "$@"
