#!/bin/bash

GENOMIC_REFERENCE_DATA=/home/s.benjamin/genomic_reference_data

main() {
  arg_parse "$@"

  SCRIPT_NAME=$(echo $0 | awk -F / '{print $NF}')
  echo \###################$SCRIPT_NAME \($(date)\)#############
  echo running: $SCRIPT_NAME "$@"
  date
  echo hostname: $(hostname)
  echo \########################################################
  echo


  select_genome

  #activate conda environment
  eval "$(micromamba shell hook --shell=bash)"
  micromamba activate /home/s.benjamin/micromamba/envs/ovation_rrbs_pipeline_2022

  time trim_illumina_adapters
  time trim_diversity_adaptors
  time align_to_genome
  ##TODO: remove PCR duplicates (optional)
  time methylation_calling
  time combine_methylation_coverage_to_tiles 100 10 #<tile_size> <min_coverage>

  echo \##############finished $SCRIPT_NAME \($(date)\)#############
}

select_genome() {
  if [[ $genome == "mm10" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/from_huji/mm10/Sequence/WholeGenomeFasta/
  elif [[ $genome == "mm9" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/
  elif [[ $genome == "hg38" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/hg38/
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
  if [[ $READ_TYPE == "single_end" ]]; then
    cmd=$(${TRIM_GALORE} --adapter AGATCGGAAGAGC $INPUT_FASTQ --cores $N_CORES --fastqc)
  else
    cmd=$(${TRIM_GALORE} --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $INPUT_FASTQ_1 $INPUT_FASTQ_2 \
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
#	echo \###################$SCRIPT_NAME \($(date)\)#############
#	echo runnig: trim_galore --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $1 $2 --cores $N_CORES --fastqc \($(date)\)
#	${TRIM_GALORE} --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC $1 $2 --cores $N_CORES --fastqc
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
#	echo \###################$SCRIPT_NAME \($(date)\)#############
#	echo runnig: trim_galore --adapter AGATCGGAAGAGC $1 --cores $N_CORES \($(date)\)
#	${TRIM_GALORE} --adapter AGATCGGAAGAGC $1  --cores $N_CORES --fastqc
#	echo \########################################################
#}
#

trim_diversity_adaptors() {
  #https://github.com/nugentechnologies/NuMetRRBS#diversity-trimming-and-filtering-with-nugens-diversity-trimming-scripts
  echo \###################$SCRIPT_NAME \($(date)\)#############
  if [[ $READ_TYPE == "single_end" ]]; then
    TRIM_GALORE_OUTPUT=$(echo $INPUT_FASTQ | awk -F / '{print $NF}' | sed 's/\(\.fastq\|.fq\)\.gz/_trimmed.fq.gz/')
    echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $TRIM_GALORE_OUTPUT
    python2 $DIVERSITY_TRIM_SCRIPT -1 $TRIM_GALORE_OUTPUT
  else
    TRIM_GALORE_OUTPUT_1=$(echo $INPUT_FASTQ_1 | awk -F / '{print $NF}' | sed 's/\(\.fastq\|.fq\)\.gz/_val_1.fq.gz/')
    TRIM_GALORE_OUTPUT_2=$(echo $INPUT_FASTQ_2 | awk -F / '{print $NF}' | sed 's/\(\.fastq\|.fq\)\.gz/_val_2.fq.gz/')
    echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $TRIM_GALORE_OUTPUT_1 -2 $TRIM_GALORE_OUTPUT_2
    python2 $DIVERSITY_TRIM_SCRIPT -1 $TRIM_GALORE_OUTPUT_1 -2 $TRIM_GALORE_OUTPUT_2
  fi

  echo \########################################################
}

align_to_genome() {
  BISMARK_GENOME_LOCATION=/utemp/s.benjamin/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta
  if [[ $READ_TYPE == "single_end" ]]; then
    TRIM_DIVERSITY_OUTPUT=$(echo $TRIM_GALORE_OUTPUT | sed 's/\.gz/_trimmed.fq.gz/')
    COMMAND=$(echo $BISMARK $non_directional --multicore $N_CORES --bowtie2 $BISMARK_GENOME_LOCATION $TRIM_DIVERSITY_OUTPUT)
  else
    TRIM_DIVERSITY_OUTPUT_1=$(echo $TRIM_GALORE_OUTPUT_1 | sed 's/\.gz/_trimmed.fq.gz/')
    TRIM_DIVERSITY_OUTPUT_2=$(echo $TRIM_GALORE_OUTPUT_2 | sed 's/\.gz/_trimmed.fq.gz/')
    COMMAND=$(echo $non_directional --multicore $N_CORES --bowtie2 $BISMARK_GENOME_LOCATION -1 $TRIM_DIVERSITY_OUTPUT_1 -2 $TRIM_DIVERSITY_OUTPUT_2)
  fi

  echo \###################$SCRIPT_NAME \($(date)\)#############
  echo $SCRIPT_NAME runnig: $COMMAND \($(date)\)
  $COMMAND
  echo \########################################################
}

# remove PCR duplicates()
# {

# }

methylation_calling() {
  if [[ $READ_TYPE == "single_end" ]]; then
    ALIGNMENT_OUTPUT=$(echo $TRIM_DIVERSITY_OUTPUT | sed 's/\.fq\.gz/_bismark_bt2.bam/')
    #By default, this mode will only consider cytosines in CpG context, but it can be extended to cytosines in any sequence context by using the option --CX
    COMMAND=$(echo bismark_methylation_extractor --multicore $N_CORES --bedGraph --buffer_size 10G --output methylation_extractor_output $ALIGNMENT_OUTPUT)
  else
    ALIGNMENT_OUTPUT=$(echo $TRIM_DIVERSITY_OUTPUT_1 | sed 's/\.fq\.gz/_bismark_bt2_pe.bam/')
    COMMAND=$(echo bismark_methylation_extractor -p --multicore $N_CORES --bedGraph --buffer_size 10G --output methylation_extractor_output $ALIGNMENT_OUTPUT)
  fi

  echo \###################$SCRIPT_NAME \($(date)\)#############
  echo $SCRIPT_NAME runnig: $COMMAND
  $COMMAND
  echo \########################################################
}

combine_methylation_coverage_to_tiles() {
  #positional args: <tile_size> <min_coverage>
  # adapted from Adam's script:
  # to make WholeGenome_100bpTiles.bed :
  # cat mm9.chrom.sizes | grep -v -P 'X|Y|M' | awk '{OFS="\t"; for (i=1; i<=$2; i+=100) print $1,i,i+100-1}' > WholeGenome_100bpTiles.bed
  TILE_SIZE=$1
  MIN_COVERAGE=$2
  cd methylation_extractor_output/
  METH_CALLING_OUTPUT=$(ls | grep cov.gz)
  # 100bp tiles variant 2: First calculate the tiles and then remove tiles with total coverage < 10
  FileOut=$(echo ${METH_CALLING_OUTPUT} | awk -v tile_size=$TILE_SIZE -F "." '{print $1 "_" tile_size "bp_tiles.bed" }')
  bedtools intersect -a /utemp/s.benjamin/genomic_reference_data/mm10_whole_genome_${TILE_SIZE}bpTiles.bed -b ${METH_CALLING_OUTPUT} -wa -wb | awk -v cov=${MIN_COVERAGE} -v tileSize=${TILE_SIZE} 'BEGIN {OFS="\t"; Prev=-1} {if ($2 == Prev) {T=T+$8+$9; M=M+$8} else {if (Prev!=-1 && T>=cov) {print PrevChr,Prev,Prev+tileSize-1,M/T};T=$8+$9; M=$8;}; Prev=$2; PrevChr=$1}' >../${FileOut}
  #nor sure what this is for
  #bedtools unionbedg -names `du -a -L | grep Tiles | awk '{print $2}' | sort | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'` -header -filler NA -i `du -a -L | grep Tiles | awk '{print $2}' | sort` > 100bpTiles_Tiles_Cov10_Tissues.bed

  echo \###################$SCRIPT_NAME \($(date)\)#############
  echo $SCRIPT_NAME runnig: combine_methylation_coverage_to_tiles $1 $2
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
      READ_TYPE="single_end"
      shift # past argument
      ;;
    -paired-end)
      READ_TYPE="paired_end"
      shift # past argument
      ;;
    -input_fastq_file)
      INPUT_FASTQ="$2"
      shift # past argument
      shift # past value
      ;;
    -paired_input_fastq_files)
      INPUT_FASTQ_1="$2"
      shift # past argument
      INPUT_FASTQ_2="$2"
      shift # past argument2
      shift # past value
      ;;
    -n_cores)
      N_CORES="$2"
      shift # past argument
      shift # past value
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
