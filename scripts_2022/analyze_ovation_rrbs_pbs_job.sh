#!/bin/bash

main()
{
	arg_parse "$@"
	set_software_paths
	
	#activate python virtual enviorment
	source ${PYTHON_ENV}/bin/activate 
		if [[ READ_TYPE == "single_end" ]]
		then
			trim_illumina_adapter_single_end $INPUT_FASTQ
		#else #if [[ READ_TYPE == "paired_end"]]
		fi
	


	deactivate #deactivate python virtual enviorment
}


set_software_paths()
{
	#cutadapt version 3.7
	#cutadapt is a python package installed in PYTHON_ENV.
	PYTHON_ENV=/home/s.benjamin/bioinformatics_software/ovation-rrbs-methyl-seq__python3-env
	
	#trim_galore version 0.6.8 (15 01 2022)
	TRIM_GALORE=/home/s.benjamin/bioinformatics_software/TrimGalore/trim_galore
	
	#Bowtie 2 version 2.2.3
	BOWTIE2=/usr/local/bin/bowtie2
	
	#samtools version 0.1.19-44428cd
	SAMTOOLS=/usr/local/samtools/bin/samtools
	
	#Bismark Version: v0.23.1
	BISMARK=/home/s.benjamin/bioinformatics_software/Bismark-0.23.1/bismark
	
	#java version "1.8.0_202"
	JAVA=/usr/local/matlab2021b/sys/java/jre/glnxa64/jre/bin/java
	
	#FastQC v0.11.9
	FASTQC=/home/s.benjamin/bioinformatics_software/FastQC/fastqc
	
	
	#add paths to executables
	ADD_TO_PATH=""
	for executable in $TRIM_GALORE $BOWTIE2 $SAMTOOLS $BISMARK $JAVA  $FASTQC
	do
		
		ADD_TO_PATH+=$(echo $executable | awk -F / 'NF{NF--};{OFS = FS; print $0}'):
	done
	echo $ADD_TO_PATH
	export PATH="$ADD_TO_PATH:$PATH"
}


trim_illumina_adapter_single_end()
{
	#first positional argument is fastq file to trim
	echo $0 runnig: ${TRIM_GALORE} --adapter AGATCGGAAGAGC $1 
	${TRIM_GALORE} --adapter AGATCGGAAGAGC $1  
}


arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
     -single-end)
        READ_TYPE="single_end"
        shift # past argument
        shift # past value
        ;;
     -paired-end)
        READ_TYPE="paired_end"
        shift # past argument
        shift # past value
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
