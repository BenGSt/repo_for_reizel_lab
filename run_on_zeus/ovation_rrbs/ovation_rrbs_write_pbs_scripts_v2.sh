#!/bin/bash

main() {
  arg_parse "$@"
  script=/home/s.benjamin/bioinformatics_software/repo_for_reizel_lab/run_on_zeus/ovation_rrbs/analyze_ovation_rrbs_v2.sh

  if [[ ! -d $OUTPUT_DIR ]]; then
    mkdir $OUTPUT_DIR
  fi
  cd $OUTPUT_DIR

  sample_list=$(ls $RAW_SAMPLES_DIR)
  for sample in $sample_list; do
    dir_name=$(echo $sample | awk -F . '{print $1}')
    if [[ $READ_TYPE == "single_end" ]]; then
      script_args=$(echo -genome $genome $non_directional -n_cores $N_CORES $READ_TYPE -input_fastq_file $RAW_SAMPLES_DIR/$sample \> $dir_name.log 2\>\&1)
    else #TODO test paired end, maybe add _1 _2 not only R1 R2
      echo sample = $sample
      r1=$(find $RAW_SAMPLES_DIR -name R1)
      r2=$(find $RAW_SAMPLES_DIR -name R2)
      script_args=$(echo -genome $genome $non_directional -n_cores $N_CORES $READ_TYPE -paired_input_fastq_files ${r1} ${r2} \> $dir_name.log 2\>\&1)
    fi

    mkdir $dir_name
    cd $dir_name
    cat <<EOF >ovation_rrbs_${dir_name}.q
#!/bin/bash
#PBS  -N  ovation_rrbs_${dir_name}
#PBS  -q  zeus_all_q
#PBS  -m  abe
#PBS  -M  s.benjamin@technion.ac.il
#PBS  -l select=1:ncpus=${N_CORES}
#PBS  -l select=mem=32gb
#PBS  -l walltime=24:00:00
PBS_O_WORKDIR=$(pwd)
cd \$PBS_O_WORKDIR


$script $script_args

EOF
    cd ..
  done
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
    -non-directional)
      non_directional="-non_directional"
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
    esac
  done
}

help() {
  cat <<EOF
  makes a directory for each sample and writes a pbs script ready to be run with qsub to it.

	-single-end
	or
	-paired-end
	
	-raw_samples_dir <path> (should contain a directory for each sample containing it's fastq files.)
	-output_dir
	-n_cores
	-genome <mm9/mm10/hg38>
	[-non-directional]

EOF
}

main "$@"
