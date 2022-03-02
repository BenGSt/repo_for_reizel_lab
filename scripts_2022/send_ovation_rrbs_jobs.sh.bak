#!/bin/bash

main()
{
	arg_parse "$@"
	cd $OUTPUT_DIR

	for sample in $(ls $RAW_SAMPLES_DIR| grep -P 'fastq|fq | grep -v md5')
	do
		dir_name=$(echo $sample | awk -F . '{print $1}')
		mkdir $dir_name
		cd $dir_name
		#TODO: send_job_from_here
		echo analyze_ovation_script_single_sample $READ_TYPE -input_fastq_file $sample | tee $dir_name.log
		cd ..
	done
}

arg_parse()
{

  while [[ $# -gt 0 ]]; do
  echo $1
    case $1 in
     -single-end)
        READ_TYPE="single-end"
        shift # past argument
        ;;
     -paired-end)
        READ_TYPE="paired-end"
        shift # past argument
        ;;
     -raw_samples_dir)
        RAW_SAMPLES_DIR="$2"
        shift # past argument
        shift # past value
        ;;
     -output_dir)
        OUTPUT_DIR="$2"
        shift # past argument
        shift # past value
        ;;		
      -h|--help)
        help
        exit 1
        ;;
    esac
  done
}


help()
{
cat << EOF
	-single-end)
	or
	-paired-end)
	-raw_samples_dir)
	-output_dir)
EOF
}


main "$@"
