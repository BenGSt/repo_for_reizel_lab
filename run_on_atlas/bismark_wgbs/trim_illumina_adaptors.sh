#!/bin/bash

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

help() {
  cat <<EOF
----------------------------------------
Project: Reizel Lab Bioinformatics Pipelines
Pipeline: Bismark WGBS
Script: trim_illumina_adaptors.sh
Author: Ben G. Steinberg
Last Update: 4 Sep 2023
----------------------------------------

Run after prepare_jobs.sh [and split_fastq.sh if fastq files are split]

USAGE: trim_illumina_adaptors.sh {-input-fastq-file <sample.fq.gz>  or  -paired-input-fastq-files <sample_R1.fq.gz> <sample_R2.fq.gz>}
                                 -output-dir <path> [-extra-trim-galore-options "multiple quoted options"]

Resources: $TRIM_JOB_CPUS cores, $TRIM_JOB_MEM RAM (defined in shared.sh)

Arguments:
-input-fastq-file <sample.fq.gz> or -paired-input-fastq-files <sample_R1.fq.gz> <sample_R2.fq.gz>
-output-dir <path> output will be written to this directory

optional:
-extra-trim-galore-options "multiple quoted options"

handy extra options from trim_galore manual
===========================================
--clip_R1 <int>         Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end
                      reads). This may be useful if the qualities were very poor, or if there is some
                      sort of unwanted bias at the 5' end. Default: OFF.

--clip_R2 <int>         Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end reads
                        only). This may be useful if the qualities were very poor, or if there is some sort
                        of unwanted bias at the 5' end. For paired-end BS-Seq, it is recommended to remove
                        the first few bp because the end-repair reaction may introduce a bias towards low
                        methylation. Please refer to the M-bias plot section in the Bismark User Guide for
                        some examples. Default: OFF.

--three_prime_clip_R1 <int>     Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end
                        reads) AFTER adapter/quality trimming has been performed. This may remove some unwanted
                        bias from the 3' end that is not directly related to adapter sequence or basecall quality.
                        Default: OFF.

--three_prime_clip_R2 <int>     Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER
                        adapter/quality trimming has been performed. This may remove some unwanted bias from
                        the 3' end that is not directly related to adapter sequence or basecall quality.
                        Default: OFF.

EOF
}

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  print_info "running: " "$script_name " "$@"
  arg_parse "$@"

  echo "If output files from previous runs exist, they will be removed as to not corrupt the current run."
  rm -fv $(find $output_dir -name "*trimmed*" -o -name "*val_[0-1]*" -o -name "*fastqc*")

  mkdir -p $output_dir
  cd $output_dir || exit 1

  if [[ $read_type == "single_end" ]]; then
    time trim_illumina_adapter_single_end $input_fastq
  else #if [[ $read_type == "paired_end" ]]
    time trim_illumina_adapter_paired_end $input_fastq_1 $input_fastq_2
  fi

  print_info "finished: " "$script_name " "$@"
}

trim_illumina_adapter_paired_end() { #<R1> <R2>
  #positional argument are  R1, R2 fastq file to trim
  #note on multicores from cutadapt manual
  #To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
  #note from trim_galore manual
  #It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
  #--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
  cmd="trim_galore --dont_gzip --paired $1 $2 --cores $TRIM_JOB_CPUS --fastqc $extra"
  echo runnig: $cmd
  $cmd
}

trim_illumina_adapter_single_end() { #<R1>
  cmd="trim_galore  $1 --dont_gzip --cores $TRIM_JOB_CPUS --fastqc $extra"
  echo runnig: $cmd
  $cmd
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
    -input-fastq-file)
      read_type="single_end"
      input_fastq="$2"
      ecgo DEBUG: input_fastq: $input_fastq
      shift # past argument
      shift # past value
      ;;
    -paired-input-fastq-files)
      read_type="paired_end"
      fastqs=$2
      input_fastq_1=${fastqs%% *} #first string
      input_fastq_2=${fastqs#* }  #second string
      echo DEBUG: input_fastq_1: $input_fastq_1
      echo DEBUG: input_fastq_2: $input_fastq_2
      shift
      shift
      ;;
    -output-dir)
      output_dir="$2"
      echo DEBUG: output_dir: $output_dir
      shift
      shift
      ;;
    -extra-trim-galore-options)
      extra="$2"
      echo DEBUG: extra: $extra
      shift
      shift
      ;;
    *)
      echo "ERROR: unknown argument \"$1\""
      help
      exit 1
      ;;
    esac
  done
}

main "$@"
