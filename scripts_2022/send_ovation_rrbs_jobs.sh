#!/bin/bash

main()
{
for x in $(ls ~/raw_sequencing_data/KKTR-FahRenegeration-rrbs/FASTQ/FGC2166/)
do
	mkdir $x
	cd $x
	#TODO: send_job_from_here
	cd ..
done
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
