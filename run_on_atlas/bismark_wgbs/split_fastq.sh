#!/usr/bin/env bash
source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  print_info "running: " "$script_name " "$@"
  arg_parse "$@"
  mkdir -p $output_dir
  cd $output_dir
  echo pwd: $PWD

  #if split directory exists, delete it
  if [[ -d split ]]; then
    rm -rfv split
    echo deleted preexisting split directory to make sure restarted jobs are clean
  fi
  mkdir split

  # Split fastq files, N_READS_PER_FILE reads per file
  if [[ $read_type == "single_end" ]]; then
    to_split=($input_fastq)
  else #if [[ $read_type == "paired_end" ]]
    to_split=($input_fastq_1 $input_fastq_2)
  fi

  n=0
  for fastq in "${to_split[@]}"; do
    mkfifo fifo$n
    pigz -p 1 -cd $fastq >fifo$n &
    echo split -dl $n_lines_per_file fifo$n split/$(echo $(basename $fastq) | sed 's/.fastq.gz\|fq.gz//')_chunk_ --additional-suffix=.fq
    split -dl $n_lines_per_file fifo$n split/$(echo $(basename $fastq) | sed 's/.fastq.gz\|fq.gz//')_chunk_ --additional-suffix=.fq &
    ((n++))
  done

  wait
  rm -v fifo*
  cd split

  for chunk in $(seq -w 00 $((n_chunks - 1))); do
    mkdir $chunk
    mv *$chunk.fq $chunk
  done

  print_info "finished: " "$script_name " "$@"
}

help() {
  cat <<EOF
----------------------------------------
Project: Reizel Lab Bioinformatics Pipelines
Pipeline: Bismark WGBS
Script: split_fastq.sh
Author: Ben G. Steinberg
Last Update: 4 Sep 2023
----------------------------------------

This script splits gzipped fastq files into <n_chunks> unzipped fastq files, each containing <n_reads_per_chunk> reads.
n_chunks, and n_reads_per_chunk must be given as arguments in advance because the script is intended to run as part of
the Bismark WGBS pipeline,on the Atlas cluster, which does not allow jobs to submit jobs. The number of chunks must be
known in advance in order to write matching submission files for splitting, trimming and aligning.

Usage:
  paired end:
    split_fastq.sh -output-dir <output_dir> -chunks <n_chunks> -reads-per-chunk <n_reads_per_chunk>
                   -paired-input-fastq-files <R1.fq.gz> <R2.fq.gz>
  single_end:
    split_fastq.sh -output-dir <output_dir> -chunks <n_chunks> -reads-per-chunk <n_reads_per_chunk>
                   -input_fastq <sample.fq.gz>

EOF
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
    -chunks)
      n_chunks=$2
      shift
      shift
      ;;
    -reads-per-chunk)
      n_lines_per_file=$((4 * $2)) #4 lines per read
      shift
      shift
      ;;
    -input-fastq-file)
      read_type="single_end"
      input_fastq="$2"
      shift # past argument
      shift # past value
      ;;
    -paired-input-fastq-files)
      read_type="paired_end"
      input_fastq_1="$2"
      shift
      input_fastq_2="$2"
      shift
      shift
      ;;
    -output-dir)
      output_dir="$2"
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
