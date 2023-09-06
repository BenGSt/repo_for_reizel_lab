#!/bin/bash

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

#TODO: to include all Cytosines use option --CX
#TODO: -keep-files #OT OB *.bedGraph.gz

help() {
  cat <<EOF
----------------------------------------
Project: Reizel Lab Bioinformatics Pipelines
Pipeline: Bismark WGBS
Script: methylation_calling.sh
Author: Ben G. Steinberg
Last Update: 4 Sep 2023
----------------------------------------

Run after deduplicate.sh, Produces coverage report, cov file (only CpGs)

USAGE: methylation_calling.sh -output-dir <path> [-extra-options "multiple double quoted options"] [-bam-dir <path>]

Resources: $METH_CALL_JOB_CPUS cores, $METH_CALL_JOB_MEM RAM (defined in shared.sh)

Arguments:
-output-dir <path> Alignment output bam file(s) are expected to be found here, output will be written to here as well.

optional:
[-extra-options "multiple double quoted options"] (e.g. --ignore <int>. See Bismark manual)
[-bam-dir <path>] Use for m-bias fix (when the bam files are in a different directory than the output directory).

EOF
}

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  arg_parse "$@"
  cd "$output_dir" || exit 1
  print_info "running: " "$script_name " "$@"

  #remove output files from previous runs
  echo "If output files from previous runs exist, they will be removed as to not corrupt the current run."
  rm -fv $(find . -name "*.cov.gz" -o -name "*.bedGraph.gz" -o -name "*OB*.txt*" -o -name "*OT*.txt*")

  call_methylation

  #rename the cov output:
  current_dir=$(pwd | awk -F'/' '{print $NF}')
  mv -v ./*.cov.gz ${current_dir}.cov.gz

  #cleanup
  rm -v $(find ./ | grep -P 'OT|OB')
  rm -v $(find . -name "*.bedGraph.gz")

  print_info "finished: " "$script_name " "$@"
}

call_methylation() {
  #if bam_dir is an empty string (i.e. not set), find the bam files in the current directory
  if [[ -z $bam_dir ]]; then
    alignment_output=$(find . -name '*bismark*deduplicated*bam')
  else
    alignment_output=$(find $bam_dir -name '*bismark*deduplicated*bam')
  fi

  #check if paired end and set flag
  echo $alignment_output | grep 'pe' >/dev/null && paired="-p" || paired=""

  #  command=$(echo bismark_methylation_extractor --bedgraph $paired $ignore_r2 --multicore $METH_CALL_INSTANCES --gzip --buffer_size $METH_CALL_BUFFER_SIZE $extra $alignment_output)
  command="bismark_methylation_extractor --bedgraph $paired --multicore $METH_CALL_INSTANCES --gzip --buffer_size $METH_CALL_BUFFER_SIZE $extra $alignment_output"
  #NOTE: samtools broken pipe and perl gzip: broken pipe errors occur. this is a knows issue and should not affect the output.
  #      an ld version of samtools (--samtools_path /Local/bfe_reizel/samtools-0.1.19/) fixes the samtools error but not the perl error, perhaps an old version of perl will work.
  #      Leaving this as is for now, the errors are not fatal and the output is fine.

  #--ample_memory speeds things up for samples over 10 million reads or so. since it may take over an hour to get going ATLAS policy holds the jobs.
  #  command=$(echo bismark_methylation_extractor --ample_memory --bedgraph $paired $ignore_r2 --multicore $METH_CALL_INSTANCES --gzip  $extra $alignment_output)
  echo
  echo $script_name runnig: $command
  echo
  $command
}

arg_parse() {
  while [[ $# -gt 0 ]]; do
    case $1 in
    -h | --help)
      help
      exit 1
      ;;
    -output-dir)
      output_dir="$2"
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
