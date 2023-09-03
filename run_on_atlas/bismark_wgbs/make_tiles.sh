#!/bin/bash
source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh

help() {
  cat <<EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh, align_to_genome.sh, methylation_calling.sh
	resources: 1 core, 30GB RAM

	-output_dir
	-genome <mm10 or hg38>
EOF
}

main() {
  # NOTE: not cleaning up tiles from previous runs, because they wii be overwritten anyway,
  # and people may want to run this script again with different parameters (e.g. tile size).
  arg_parse "$@"
  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
  cd $output_dir || exit 1
  script_name=$(echo $0 | awk -F / '{print $NF}')

  print_info "running: " "$script_name " "$@"
  time combine_methylation_coverage_to_tiles 100 10 $genome #<tile_size> <min_coverage> <genome>
  print_info "finished: " "$script_name " "$@"
}

combine_methylation_coverage_to_tiles() {
  #positional args: <tile_size> <min_coverage> <genome = mm10 or hg38>
  # adapted from Adam's script:
  # to make WholeGenome_100bpTiles.bed :
  # cat mm9.chrom.sizes | grep -v -P 'X|Y|M' | awk '{OFS="\t"; for (i=1; i<=$2; i+=100) print $1,i,i+100-1}' > WholeGenome_100bpTiles.bed
  # cat hg38.chrom.sizes  | awk '{OFS="\t"; for (i=1; i<=$2; i+=100) print $1,i,i+100-1}' > hg38_100bp_tiles.bed

  tile_size=$1
  min_coverage=$2
  #	cd methylation_extractor_output/ #trying without this folder
  meth_calling_output=$(find . -name "*.cov.gz")

  mm10_tiles=${GENOMIC_REFERENCE_LOCATION}/mm10_whole_genome_${tile_size}bpTiles.bed
  hg38_tiles=${GENOMIC_REFERENCE_LOCATION}/hg38/hg38_100bp_tiles.bed #TODO: create tiles for minChromSet, speed things up.
  if [[ $3 == "mm10" ]]; then
    tiles_file=$mm10_tiles
  elif [[ $3 == "hg38" ]]; then
    tiles_file=$hg38_tiles
  fi

  # 100bp tiles variant 2: First calculate the tiles and then remove tiles with total coverage < 10
  #	output_file=$(echo ${meth_calling_output} | awk -v tile_size=$tile_size -F "." '{print $1 "_" tile_size "bp_tiles.bed" }')
  #	output_file=$(echo ${meth_calling_output} | sed 's/_R[1-2].*/_100bp_tiles.bed/')

  #get current directory name:
  current_dir=$(pwd | awk -F'/' '{print $NF}')
  output_file=${current_dir}_100bp_tiles.bed
  bedtools intersect -a $tiles_file -b $meth_calling_output -wa -wb | awk -v cov=$min_coverage -v tileSize=$tile_size 'BEGIN {OFS="\t"; Prev=-1} {if ($2 == Prev) {T=T+$8+$9; M=M+$8} else {if (Prev!=-1 && T>=cov) {print PrevChr,Prev,Prev+tileSize-1,M/T};T=$8+$9; M=$8;}; Prev=$2; PrevChr=$1}' >${output_file}

  #to unite all tiles from different samples:
  #bedtools unionbedg -names `du -a -L | grep Tiles | awk '{print $2}' | sort | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'` -header -filler NA -i `du -a -L | grep Tiles | awk '{print $2}' | sort` > 100bpTiles_Tiles_Cov10_Tissues.bed
}

arg_parse() {
  while [[ $# -gt 0 ]]; do
    case $1 in
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
    -* | --*)
      help
      exit 1
      ;;
    -h | --help)
      help
      exit 1
      ;;
    esac
  done
}

main "$@"
