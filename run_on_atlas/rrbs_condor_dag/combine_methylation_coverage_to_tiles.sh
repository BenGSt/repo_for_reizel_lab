#!/bin/bash

GENOMIC_REFERENCE_LOCATION=/storage/bfe_reizel/bengst/genomic_reference_data
BISMARK_GENOME_LOCATION=${GENOMIC_REFERENCE_LOCATION}/from_huji/mm10/Sequence/WholeGenomeFasta

help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh, align_to_genome.sh, methylation_calling.sh
	resources: 1 core, 3GB RAM

	-input_fastq_file <sample.fq.gz> or -paired_input_fastq_files <sample_R1.fq.gz> <sample_R2.fq.gz>
	-output_dir
EOF
}


main()
{
  arg_parse "$@"
	set_software_paths
	cd $output_dir
	script_name=$(echo $0 | awk -F / '{print $NF}')

	echo
	echo
	echo \################################
	echo \################################
	echo running: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo

  #TODO: rm bam file and only keep bismark.cov?
	time combine_methylation_coverage_to_tiles 100 10 #<tile_size> <min_coverage>

	echo
	echo
	echo \################################
	echo \################################
	echo finished: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo
}


combine_methylation_coverage_to_tiles()
{
	#positional args: <tile_size> <min_coverage>
	# adapted from Adam's script:
	# to make WholeGenome_100bpTiles.bed :
	# cat mm9.chrom.sizes | grep -v -P 'X|Y|M' | awk '{OFS="\t"; for (i=1; i<=$2; i+=100) print $1,i,i+100-1}' > WholeGenome_100bpTiles.bed
	
	echo $script_name runnig: combine_methylation_coverage_to_tiles $1 $2
	
	tile_size=$1
	min_coverage=$2
	cd methylation_extractor_output/
	meth_calling_output=$(ls |grep cov.gz)

	# 100bp tiles variant 2: First calculate the tiles and then remove tiles with total coverage < 10
	output_file=$(echo ${meth_calling_output} | awk -v tile_size=$tile_size -F "." '{print $1 "_" tile_size "bp_tiles.bed" }')
	bedtools intersect -a ${GENOMIC_REFERENCE_LOCATION}/mm10_whole_genome_${tile_size}bpTiles.bed -b ${meth_calling_output} -wa -wb | awk -v cov=${min_coverage} -v tileSize=${tile_size} 'BEGIN {OFS="\t"; Prev=-1} {if ($2 == Prev) {T=T+$8+$9; M=M+$8} else {if (Prev!=-1 && T>=cov) {print PrevChr,Prev,Prev+tileSize-1,M/T};T=$8+$9; M=$8;}; Prev=$2; PrevChr=$1}' > ../${output_file}

	#to unite all tiles from different samples:
	#bedtools unionbedg -names `du -a -L | grep Tiles | awk '{print $2}' | sort | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'` -header -filler NA -i `du -a -L | grep Tiles | awk '{print $2}' | sort` > 100bpTiles_Tiles_Cov10_Tissues.bed
}


set_software_paths()
{
  #	echo setting software paths
	# on Atlas I'm using conda for env managment, so the software should be in path already.
	# leaving this function if I want to change something without the env.
	#TODO: consider using calls to program name from the rest of the script and delete this function.

  source /Local/bfe_reizel/anaconda3/bin/activate ovation_rrbs_pipeline_2022

	PYTHON3=python

	TRIM_GALORE=trim_galore

	BOWTIE2=bowtie2

	SAMTOOLS=samtools

	BISMARK=bismark

	FASTQC=fastqc

	DIVERSITY_TRIM_SCRIPT=/Local/bfe_reizel/anaconda3/envs/ovation_rrbs_pipeline_2022/bin/trimRRBSdiversityAdaptCustomers.py

	BEDTOOLS=bedtools

	#tecan nudup tool for pcr duplicates
	NUDUP=/Local/bfe_reizel/anaconda3/envs/ovation_rrbs_pipeline_2022/bin/nudup.py
}




arg_parse()
{
  #Note: don't need any pf these except for output_dir, the rest are left here for compatibility to rrbs_jobs.args
  #TODO: adjust this so it stayes commpatible with rrbs_jobs.args but doesn't parse unnecessary args.
  while [[ $# -gt 0 ]]; do
    case $1 in
     -single-end)
        read_type="single_end"
        shift # past argument
        ;;
     -paired-end)
        read_type="paired_end"
        shift # past argument
        ;;
     -input_fastq_file)
        input_fastq="$2"
        shift # past argument
        shift # past value
        ;;
      -paired_input_fastq_files)
        input_fastq_1="$2"
        shift # past argument
        input_fastq_2="$2"
        shift # past argument2
        shift # past value
        ;;
	-output_dir)
        output_dir="$2"
        shift # past argument
        shift # past value
        ;;
      -*|--*)
        help
        exit 1
        ;;
      -h|--help)
        help
        exit 1
        ;;
    esac
  done
}


main "$@"