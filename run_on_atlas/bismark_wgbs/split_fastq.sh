#!/usr/bin/env bash

main() {
  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
  arg_parse "$@"
  mkdir -p $output_dir
  cd $output_dir
  mkdir split

  # Split fastq files, N_READS_PER_FILE reads per file
  if [[ $read_type == "single_end" ]]; then
    to_split=($input_fastq)
  else #if [[ $read_type == "paired_end" ]]
    to_split=($input_fastq_1 $input_fastq_2)
  fi

  for fastq in "${to_split[@]}"; do
    split -dl $n_lines_per_file <(pigz -p 4 -cd $fastq) split/$(echo $(basename $fastq) | sed 's/.fastq.gz\|fq.gz//')_chunk_ --additional-suffix=.fq
  done

  cd split
  for chunk in $(seq -w 00 $n_chunks); do
    mkdir $chunk
    mv *$chunk.fq $chunk
  done
}

help() {
  #TODO: add help
  cat <<EOF
  Usage: #TODO: add usage
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