#!/usr/bin/env bash
# USAGE: fix_mbias.sh --biased_dir <path> --ignore_r1 <int> --ignore_r2  --ignore_3prime <int> --ignore_3prime_r2 <int> [--output-dir <output_dir>]

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

save_cmd() {
  if [[ $# -gt 2 ]]; then #don't (re)write cmd.txt if no args
    echo This command was run to fix m-bias. the original command was: >cmd.txt
    printf "# "
    cat $biased_dir/cmd.txt >>cmd.txt
    echo >>cmd.txt
    echo >>cmd.txt
    echo the m-bias fix command was: >>cmd.txt
    echo "$0" "$@" >>cmd.txt #TODO: preserve quotes that may be in args
  fi
}

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  print_info "running: " "$script_name " "$@"

  bias_fix=1 #for write_sample_dag_file()
  arg_parse "$@"
  extra_meth_opts="-extra-options '$ignore_r1 $ignore_r2 $ignore_3prime $ignore_3prime_r2'"
  if [[ ! $output_dir ]]; then
    echo output_dir=${biased_dir}_mbias_fixed
    output_dir=${biased_dir}_mbias_fixed
  fi
  mkdir -p $output_dir
  cd $output_dir
  save_cmd "$@"
  mkdir -p logs

  #write bismark_methylation_extractor sub files
  #find the sample directories (the first node in the path) for which bam files exist
  for sample_name in $(find $biased_dir -name "*deduplicated*bam" | awk -F / '{print $(NF-1)}'); do
    {
      unset split sep chunk
      sample_names+=($sample_name)
      mkdir -p $sample_name
      mkdir -p condor_submission_files/$sample_name
      mkdir -p logs/$sample_name

      # the following sub files are not dependent on splitting
      write_methylation_calling_job_submission_file "-bam-dir $biased_dir/$sample_name"
      write_bam2nuc_job_submission_file -override_genome $(find $biased_dir/ -name "*make_tiles.out" | head -1 | xargs grep -o -- "-genome.*" | awk '{print $2}')
      write_make_tiles_job_submission_file #genome already set by previous -override_genome
      write_multiqc_job_submission_file
      write_sample_dag_file
    }
  done

  write_top_level_dag -mbias-fix

  #list jobs and the commands to run them
  #ask if user wants to run them now, if so, run them.
  echo Submit all jobs by running: condor_submit_dag $output_dir/condor_submission_files/submit_all_bismark_wgbs.dag
  echo
  echo To run samples individually:
  for dag in $(find $output_dir/condor_submission_files -name "*.dag" | grep -v submit_all); do
    echo condor_submit_dag $dag
  done
  echo
  printf 'Submit all jobs now? (y/n) '
  read answer
  if [ "$answer" != "${answer#[Yy]}" ]; then # this grammar (the #[] operator) means that the variable $answer where any Y or y in 1st position will be dropped if they exist.
    condor_submit_dag $output_dir/condor_submission_files/submit_all_bismark_wgbs.dag
  fi
  echo Good Luck!
  cat <<EOF

Unless you need them, it is recommended to delete the bam files when you are done.

To do so, run: rm -v $(find $biased_dir -name "*deduplicated*bam")

Please download your data and delete it from atlas as soon as you are done.
!Good luck and happy clustering!
EOF

    print_info "finished: " "$script_name " "$@"
}

help() {
  cat <<EOF
  $script_name --biased_dir <path> {at least one of: --ignore_r1 <int> --ignore_r2  --ignore_3prime <int> --ignore_3prime_r2 <int>} [--output-dir <output_dir>]
obligatory options:
   --biased_dir <path>
at least one of the following:
   --ignore <int>
   --ignore_3prime <int>

   --ignore_r2 <int>
   --ignore_3prime_r2 <int>
non-obligatory options:
   [--output-dir <output_dir>] defaults to \${biased_dir}_mbias_fixed

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
      output_dir=$(realpath $2)
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
