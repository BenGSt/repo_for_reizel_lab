#!/bin/bash

main()
{
	arg_parse "$@"
	
	script=/home/s.benjamin/bioinformatics_software/repo_for_reizel_lab/scripts_2022/analyze_ovation_rrbs.sh
	
	if [[ ! -d $OUTPUT_DIR ]]; then
		mkdir $OUTPUT_DIR
	fi
	
	cd $OUTPUT_DIR

	for sample in $(ls $RAW_SAMPLES_DIR| grep -P 'fastq|fq' | grep -v md5)
	do
		dir_name=$(echo $sample | awk -F . '{print $1}')
		mkdir $dir_name
		cd $dir_name
		
		cat << EOF > ovation_rrbs_${dir_name}.q
#!/bin/sh
#PBS  -N  ovation_rrbs_${dir_name}
#PBS  -q  zeus_all_q 
#PBS  -m  abe
#PBS  -M  s.benjamin@technion.ac.il
#PBS  -l select=1:ncpus=${N_CORES}
#PBS  -l select=mem=64gb
#PBS  -l walltime=24:00:00
PBS_O_WORKDIR=$dir_name
cd \$PBS_O_WORKDIR


$script -n_cores $N_CORES $READ_TYPE -input_fastq_file $RAW_SAMPLES_DIR/$sample >  $dir_name.log 2>&1
EOF
		cd ..
	done
}

arg_parse()
{

  while [[ $# -gt 0 ]]; do
    case $1 in
     -single-end)
        READ_TYPE="-single-end"
        shift # past argument
        ;;
     -paired-end)
        READ_TYPE="-paired-end"
        shift # past argument
        ;;
     -raw_samples_dir)
        RAW_SAMPLES_DIR="$(realpath $2)"
        shift # past argument
        shift # past value
        ;;
     -output_dir)
        OUTPUT_DIR="$2"
        shift # past argument
        shift # past value
        ;;		
	 -n_cores)
        N_CORES="$2"
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
	-n_cores)
EOF
}


main "$@"
