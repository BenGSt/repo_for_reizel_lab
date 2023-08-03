#!/bin/bash

#NOTE: expects to be run in the dmrs_pipeline_2022 conda environment
GENOMIC_REFERENCE_DATA=/home/s.benjamin/genomic_reference_data

#DEFAULTS:
n_cores=10
reference_genome=/storage/bfe_reizel/bengst/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta/genome.fa
output_dir=$PWD

main() {
  arg_parse "$@"
  mkdir -p $output_dir || exit 1
  cd $output_dir || exit 1
  select_genome

  if [[ $bg_bed_file == "" ]]; then
    shuf -n 50000 $make_bg_from >random_50000_100bp_tiles_for_homer_bg.bed
    bg_bed_file=random_50000_100bp_tiles_for_homer_bg.bed
  fi

  fasta_bg=$(basename $(echo $bg_bed_file | sed 's/.bed/.fa/'))
  results_dir=$(basename $(echo $regions_bed_file | sed s/".bed"/"_HOMER_RESULTS"/))
  regions_fasta=$(basename $(echo $regions_bed_file | sed s/".bed"/".fa"/))

  printf "\n\n";
  echo fasta_bg: $fasta_bg
  echo results_dir: $results_dir
  echo regions_fasta: $regions_fasta
  printf "\n\n";

  #get bg fasta
  bedtools getfasta -fi $reference_genome -bed $bg_bed_file -fo $fasta_bg
  #e.g.: bedtools getfasta -fi $reference_genome -bed ./dmrs_plus_random_50000_100bp_tiles.bed  -fo dmrs_plus_random_50000_100bp_tiles.fa

  #Extract Sequence from bed file:
  bedtools getfasta -fi $reference_genome -bed $regions_bed_file -fo $regions_fasta

  #find motifs with  background
  findMotifs.pl $regions_fasta fasta $results_dir -p $n_cores -fastaBg $fasta_bg
}

select_genome() {
  if [[ $genome == "mm10" ]]; then
    reference_genome=$GENOMIC_REFERENCE_DATA/from_huji/mm10/Sequence/WholeGenomeFasta/genome.fa
  elif [[ $genome == "mm9" ]]; then
    reference_genome=$GENOMIC_REFERENCE_DATA/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
  elif [[ $genome == "hg38" ]]; then
    reference_genome=$GENOMIC_REFERENCE_DATA/hg38/hg38.analysisSet.fa
  else
    echo ERROR: genome $genome not supported
    exit 1
  fi
}

help() {
  cat <<EOF
#NOTE: expects to be run in the dmrs_pipeline_2022 conda environment

USAGE: $0 -regions_bed_file <regions_file.bed> {-make_bg_from <all_100bp.bed>  or -bg_bed_file <bg_file.bed>} [-genome <mm10, mm9, hg38>] [-n_cores <10>]  [-output_dir <path>]

-regions_bed_file: bed file with regions to find motifs in
-genome: mm10, mm9, hg38 (default: mm10)
-n_cores: number of cores to use (default: 10)
-make_bg_from <all_100bp.bed>: if given, make bg file from given bed file (dmr pipeline: all_100bp_tiles_united.bed

[-bg_bed_file: if -make_bg_from is not used, a bed file with background regions may be provided (relevant random regions, e.g. for DMRs based on RRBS use the DMRs + random tiles coverd by rrbs. If using as part of the dmrs_pipeline use the random_50000_100bp_tiles_for_homer_bg.bed file]
EOF
}

arg_parse() {
  [[ $# -eq 0 ]] && {
    help
    exit 1
  }
  while [[ $# -gt 0 ]]; do
    case $1 in
    -output_dir)
      output_dir=$(realpath $2)
      shift
      shift
      ;;
    -make_bg_from)
      make_bg_from=$(realpath $2)
      shift
      shift
      ;;
    -regions_bed_file)
      regions_bed_file=$(realpath $2)
      shift
      shift
      ;;
    -bg_bed_file)
      bg_bed_file=$(realpath $2)
      shift
      shift
      ;;
    -n_cores)
      n_cores="$2"
      shift
      shift
      ;;
    -genome)
      genome=$2
      shift
      shift
      ;;
    -h | --help)
      help
      exit 1
      ;;
    -GENOMIC_REFERENCE_DATA)
      GENOMIC_REFERENCE_DATA=$2
      shift
      shift
      ;;
    esac

  done
}

main "$@"
