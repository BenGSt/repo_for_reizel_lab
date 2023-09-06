#!/bin/bash

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

help() {
  cat <<EOF
----------------------------------------
Project: Reizel Lab Bioinformatics Pipelines
Pipeline: Bismark WGBS
Script: nucleotide_coverage_report.sh
Author: Ben G. Steinberg
Last Update: 4 Sep 2023
----------------------------------------

Run after methylation_calling.sh, Checks for bias by comparing the mono- and di-nucleotide coverage of the reads to the genome,
and writes an HTML report regarding the full bismark suite.

USAGE: nucleotide_coverage_report.sh -output-dir <path> -genome <hg38 or mm10>

ARGUMENTS:
  -output-dir <sample_dir>
  -genome <hg38 or mm10>

From Bismark Docs:
===================
The script bam2nuc reads BAM files and calculates the mono- and di-nucleotide coverage of the reads
(using the genomic sequence rather than the observed sequence in the reads themselves)
and compares it to the average genomic sequence composition. Reads harbouring InDels are not taken into consideration.
Mono- or dinucleotides containing Ns are ignored as well.

EOF
}

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  print_info "running: " "$script_name " "$@"

  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
  arg_parse "$@"
  cd "$output_dir" || exit 1

  if [[ $genome == "mm10" ]]; then
    bismark_genome_location=$MM10_REF # (defined in shared.sh)
  elif [[ $genome == "hg38" ]]; then
    bismark_genome_location=$HG38_REF
  else
    echo genome not recognized
    exit 1
  fi

  #remove output files from previous runs
  echo "If output files from previous runs exist, they will be removed as to not corrupt the current run."
  rm -fv $(find . -name "nucleotide_stats.txt" -o -name "report.html")

  bam2nuc --genome_folder $bismark_genome_location ./*.bam
  bismark2report --splitting_report *splitting_report.txt --mbias_report *M-bias.txt \
    --nucleotide_report *nucleotide_stats.txt --dedup_report *deduplication_report.txt

  print_info "finished: " "$script_name " "$@"
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
    -genome)
      genome=$2
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
