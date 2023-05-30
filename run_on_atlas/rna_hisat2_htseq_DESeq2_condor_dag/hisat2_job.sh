#!/bin/bash
#r1=$1
#r2=$2
#output_file=$3
#summary_file=$4 # {}_hisat2.summary.txt

main()
{
  arg_parse "$@"

	if [[ $genome == "mm10" ]]; then
	 hisat2_idx=/storage/bfe_reizel/bengst/genomic_reference_data/from_huji/mm10/Sequence/whole_genome_fasta_hisat2_index/hisat2-index
	elif [[ $genome == "hg38" ]]; then
	  hisat2_idx=/storage/bfe_reizel/bengst/genomic_reference_data/hg38/analysisSet/chromosomes_hisat2_index/hisat2-index
	else
	  echo genome not recognized
	  exit 1
  fi

  source /Local/bfe_reizel/anaconda3/bin/activate rna-seq_hisat2_htseq_deseq2_2022
  mkdir -p $(dirname $summary_file)

  if [[ $read_type == "single_end" ]]; then
    hisat2 -p 10 -x $hisat2_idx -U $r1 --summary-file $summary_file | samtools sort -n --output-fmt BAM > $output_file
  else
    hisat2 -p 10 -x $hisat2_idx -1 $r1 -2 $r2 --summary-file $summary_file | samtools sort -n --output-fmt BAM > $output_file
  fi
}
    #hisat2 -p 10 -x $hisat2_idx -1 $r1 -2 $r2 --summary-file $summary_file | samtools sort -n --output-fmt BAM >$output_file






#TODO: document options

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
      read_type="single_end"
      shift
      ;;
    -paired-end)
      read_type="paired_end"
      shift
      ;;
    -genome)
      genome="$2"
      shift
      shift
      ;;
    -r1)
      r1="$2"
      shift
      shift
      ;;
    -r2)
      r2="$2"
      shift
      shift
      ;;
    -output-file)
      output_file="$2"
      shift
      shift
      ;;
    -summary-file)
      summary_file="$2" # {}_hisat2.summary.txt
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