#!/bin/bash

GENOMIC_REFERENCE_DATA=/home/s.benjamin/genomic_reference_data
DIVERSITY_TRIM_SCRIPT=/home/s.benjamin/bioinformatics_software/NuMetRRBS/trimRRBSdiversityAdaptCustomers.py
N_CORES=20
TILE_SIZE=100
TILE_MIN_COV=10

main() {
  #enable "exit immediately on error" behavior
  set -e
  arg_parse "$@"

  #activate conda environment
  eval "$(micromamba shell hook --shell=bash)"
  micromamba activate /home/s.benjamin/micromamba/envs/ovation_rrbs_pipeline_2022

  script_name=$(echo $0 | awk -F / '{print $NF}')
  echo \###################$script_name \($(date)\)#############
  echo running: $script_name "$@"
  date
  echo hostname: $(hostname)
  echo \########################################################
  echo

  select_genome
  trim_illumina_adapters
  trim_diversity_adaptors
  align_to_genome
  ##TODO: remove PCR duplicates (optional)
  methylation_calling
  combine_methylation_coverage_to_tiles $TILE_SIZE $TILE_MIN_COV

  echo \########################################################
  echo finished $script_name \($(date)\)
}

select_genome() {
  if [[ $genome == "mm10" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/from_huji/mm10/Sequence/WholeGenomeFasta/
    genome_tiles=/home/s.benjamin/genomic_reference_data/from_huji/mm10/mm10_${TILE_SIZE}bp_tiles.bed
  elif [[ $genome == "mm9" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/
    genome_tiles=/home/s.benjamin/genomic_reference_data/Mus_musculus/UCSC/mm9/mm9_${TILE_SIZE}bp_tiles.bed
  elif [[ $genome == "hg38" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/hg38/
    genome_tiles=/home/s.benjamin/genomic_reference_data/hg38/hg38_${TILE_SIZE}bp_tiles.bed
  else
    echo ERROR: genome $genome not supported
    exit 1
  fi
}

#:'
#Positional argument are  R1, R2 fastq file to trim.
#Note on multicore from cutadapt manual:
#  To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
#Note from trim_galore manual:
#  It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
#  --cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
#'
trim_illumina_adapters() {
  if [[ $read_type == "single_end" ]]; then
    cmd=$(echo trim_galore --adapter AGATCGGAAGAGC $input_fastq --cores $N_CORES --fastqc)
  else
    cmd=$(echo trim_galore --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $input_fastq_1 $input_fastq_2 \
      --cores $N_CORES --fastqc)
  fi
  echo runnig: $cmd
  $cmd
}

#TODO: test new trim function, if works delete next two functions.
#trim_illumina_adapter_paired_end() #<R1> <R2>
#{
#	#positional argument are  R1, R2 fastq file to trim
#	#note on nulticores from cutadapt manual
#		#To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
#	#note from trim_galore manual
#		#It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
#		#--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
#	echo \###################$script_name \($(date)\)#############
#	echo runnig: trim_galore --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $1 $2 --cores $N_CORES --fastqc \($(date)\)
#	trim_galore --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $1 $2 --cores $N_CORES --fastqc
#	echo \########################################################
#}
#
#trim_illumina_adapter_single_end()
#{
#	#first positional argument is fastq file to trim
#	#note on nulticores from cutadapt manual
#		#To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
#	#note from trim_galore manual
#		#It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
#		#--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
#	echo \###################$script_name \($(date)\)#############
#	echo runnig: trim_galore --adapter AGATCGGAAGAGC $1 --cores $N_CORES \($(date)\)
#	trim_galore --adapter AGATCGGAAGAGC $1  --cores $N_CORES --fastqc
#	echo \########################################################
#}
#

trim_diversity_adaptors() {
  #https://github.com/nugentechnologies/NuMetRRBS#diversity-trimming-and-filtering-with-nugens-diversity-trimming-scripts
  if [[ $read_type == "single_end" ]]; then
    trim_galore_output=$(echo $input_fastq | awk -F / '{print $NF}' | sed 's/\(\.fastq\|.fq\)\.gz/_trimmed.fq.gz/')
    echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output
    python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output
  else
    trim_galore_output_1=$(echo $input_fastq_1 | awk -F / '{print $NF}' | sed 's/\(\.fastq\|.fq\)\.gz/_val_1.fq.gz/')
    trim_galore_output_2=$(echo $input_fastq_2 | awk -F / '{print $NF}' | sed 's/\(\.fastq\|.fq\)\.gz/_val_2.fq.gz/')
    echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output_1 -2 $trim_galore_output_2
    python2 $DIVERSITY_TRIM_SCRIPT -1 $trim_galore_output_1 -2 $trim_galore_output_2
  fi
}

align_to_genome() {
  if [[ $read_type == "single_end" ]]; then
    trim_diversity_output=$(echo $trim_galore_output | sed 's/\.gz/_trimmed.fq.gz/')
    cmd=$(echo bismark $non_directional --parallel $(($N_CORES / 4)) --bowtie2 $bismark_genome_location $trim_diversity_output)
  else
    trim_diversity_output_1=$(echo $trim_galore_output_1 | sed 's/\.gz/_trimmed.fq.gz/')
    trim_diversity_output_2=$(echo $trim_galore_output_2 | sed 's/\.gz/_trimmed.fq.gz/')
    cmd=$(echo $non_directional --parallel $(($N_CORES / 4)) --bowtie2 $bismark_genome_location -1 $trim_diversity_output_1 -2 $trim_diversity_output_2)
  fi

  echo runnig: $cmd \($(date)\)
  $cmd

}

# remove PCR duplicates()
# {

# }

methylation_calling() {
  if [[ $read_type == "single_end" ]]; then
    alignment_output=$(echo $trim_diversity_output | sed 's/\.fq\.gz/_bismark_bt2.bam/')
    #By default, this mode will only consider cytosines in CpG context, but it can be extended to cytosines in any sequence context by using the option --CX
    cmd=$(echo bismark_methylation_extractor --multicore $(($N_CORES / 3)) --bedGraph --buffer_size 10G --output methylation_extractor_output $alignment_output)
  else
    alignment_output=$(echo $trim_diversity_output_1 | sed 's/\.fq\.gz/_bismark_bt2_pe.bam/')
    cmd=$(echo bismark_methylation_extractor -p --multicore $(($N_CORES / 3)) --bedGraph --buffer_size 10G --output methylation_extractor_output $alignment_output)
  fi

  echo runnig: $cmd
  $cmd

}

combine_methylation_coverage_to_tiles() {
  #positional args: <tile_size> <min_coverage>
  # adapted from Adam's script:
  # to make WholeGenome_100bpTiles.bed :
  # cat mm9.chrom.sizes | grep -v -P 'X|Y|M' | awk '{OFS="\t"; for (i=1; i<=$2; i+=100) print $1,i,i+100-1}' > WholeGenome_100bpTiles.bed
  tile_size=$1
  min_coverage=$2

  cd methylation_extractor_output/
  meth_call_output=$(ls | grep cov.gz)

  # 100bp tiles variant 2: First calculate the tiles and then remove tiles with total coverage < 10
  file_out=$(echo meth_call_output | awk -v tile_size=$tile_size -F "." '{print $1 "_" tile_size "bp_tiles.bed" }')
  bedtools intersect -a $genome_tiles -b $meth_call_output -wa -wb | awk -v cov=$min_coverage -v tileSize=$tile_size 'BEGIN {OFS="\t"; Prev=-1} {if ($2 == Prev) {T=T+$8+$9; M=M+$8} else {if (Prev!=-1 && T>=cov) {print PrevChr,Prev,Prev+tileSize-1,M/T};T=$8+$9; M=$8;}; Prev=$2; PrevChr=$1}' >../${file_out}

  #unite tiles from different samples into big table
  #bedtools unionbedg -names `du -a -L | grep Tiles | awk '{print $2}' | sort | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'` -header -filler NA -i `du -a -L | grep Tiles | awk '{print $2}' | sort` > 100bpTiles_Tiles_Cov10_Tissues.bed

  echo \###################$script_name \($(date)\)#############
  echo $script_name runnig: combine_methylation_coverage_to_tiles $1 $2
  echo \########################################################

}

help() {
  cat <<EOF
	-single-end or -paired-end
	-input_fastq_file or -paired_input_fastq_files
	-n_cores
EOF
}

arg_parse() {
  while [[ $# -gt 0 ]]; do
    case $1 in
    -single-end)
      read_type="single_end"
      shift
      ;;
    -paired-end)
      read_type="paired_end"
      shift
      ;;
    -input_fastq_file)
      input_fastq="$2"
      shift # past argument
      shift # past value
      ;;
    -paired_input_fastq_files)
      input_fastq_1="$2"
      shift
      input_fastq_2="$2"
      shift 2
      shift
      ;;
    -n_cores)
      N_CORES="$2"
      shift
      shift
      ;;
    -non-directional)
      non_directional="--non_directional"
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
