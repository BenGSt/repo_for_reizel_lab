#!/bin/bash

main()
{
	SCRIPT_NAME=$(echo $0 | awk -F / '{print $NF}')
	echo \###################$SCRIPT_NAME \($(date)\)#############
	echo running: $SCRIPT_NAME "$@"
	date
	echo hostname: $(hostname)
	echo \########################################################
	echo
	
	arg_parse "$@"
	set_software_paths
	
	#activate python virtual enviorment
	source ${PYTHON_ENV}/bin/activate 
		if [[ $READ_TYPE == "single_end" ]]
		then
			time trim_illumina_adapter_single_end $INPUT_FASTQ
			time trim_diversity_adaptors
			time align_to_genome
			
		#else #if [[ READ_TYPE == "paired_end" ]]
			#TODO
		fi
	deactivate #deactivate python virtual enviorment
	echo \##############finished $SCRIPT_NAME \($(date)\)#############
}


trim_illumina_adapter_single_end()
{
	#first positional argument is fastq file to trim
	#note on nulticores from cutadapt manual
		#To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
	#note from trim_galore manual 
		#It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
		#--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
	echo \###################$SCRIPT_NAME \($(date)\)#############
	echo runnig: trim_galore --adapter AGATCGGAAGAGC $1 --cores 4 \($(date)\)
	${TRIM_GALORE} --adapter AGATCGGAAGAGC $1  --cores 4 --fastqc
	echo \########################################################
}


trim_diversity_adaptors()
{
	#https://github.com/nugentechnologies/NuMetRRBS#diversity-trimming-and-filtering-with-nugens-diversity-trimming-scripts
	echo \###################$SCRIPT_NAME \($(date)\)#############
	TRIM_GALORE_OUTPUT=$(echo $INPUT_FASTQ |awk -F / '{print $NF}'| sed 's/\.fastq\.gz/_trimmed.fq.gz/')
	echo runnig: python2 $DIVERSITY_TRIM_SCRIPT -1 $TRIM_GALORE_OUTPUT 
	
	python2 $DIVERSITY_TRIM_SCRIPT -1 $TRIM_GALORE_OUTPUT
	echo \########################################################
}


align_to_genome()
{
	echo \###################$SCRIPT_NAME \($(date)\)#############
	TRIM_DIVERSITY_OUTPUT=$(echo $TRIM_GALORE_OUTPUT | sed 's/\.gz/_trimmed.fq.gz/')
	BISMARK_GENOME_LOCATION=/home/s.benjamin/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta
	echo $SCRIPT_NAME runnig: bismark --multicore 4 --bowtie2 $BISMARK_GENOME_LOCATION $TRIM_DIVERSITY_OUTPUT \($(date)\)
	${BISMARK} --multicore 4 --bowtie2 $BISMARK_GENOME_LOCATION $TRIM_DIVERSITY_OUTPUT
	echo \########################################################
}

set_software_paths()
{

	echo $SCRIPT_NAME: setting path \($(date)\)
	#cutadapt version 3.7
	#cutadapt is a python package installed in PYTHON_ENV.
	PYTHON_ENV=/home/s.benjamin/bioinformatics_software/ovation-rrbs-methyl-seq__python3-env
	
	#java version "1.8.0_202"
	JAVA=/usr/local/matlab2021b/sys/java/jre/glnxa64/jre/bin/java
	
	#pigz (multicore gzip needed for cutadapt multicore support)
	PIGZ=/home/s.benjamin/other_software/pigz/pigz
	

	#trim_galore version 0.6.8 (15 01 2022)
	TRIM_GALORE=/home/s.benjamin/bioinformatics_software/TrimGalore/trim_galore
	
	#Bowtie 2 version 2.2.3
	BOWTIE2=/usr/local/bin/bowtie2
	
	#samtools version 0.1.19-44428cd
	SAMTOOLS=/usr/local/samtools/bin/samtools
	
	#Bismark Version: v0.23.1
	BISMARK=/home/s.benjamin/bioinformatics_software/Bismark-0.23.1/bismark
	
	#FastQC v0.11.9
	FASTQC=/home/s.benjamin/bioinformatics_software/FastQC/fastqc
	
	#nugen diversity trimming script
	DIVERSITY_TRIM_SCRIPT=/home/s.benjamin/bioinformatics_software/NuMetRRBS/trimRRBSdiversityAdaptCustomers.py
	
	
	#add paths to executables
	ADD_TO_PATH=""
	for executable in $JAVA $PIGZ $TRIM_GALORE $BOWTIE2 $SAMTOOLS $BISMARK $FASTQC
	do
		
		ADD_TO_PATH+=$(echo $executable | awk -F / 'NF{NF--};{OFS = FS; print $0}'):
	done
	export PATH="$ADD_TO_PATH:$PATH"
}


arg_parse()
{
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
