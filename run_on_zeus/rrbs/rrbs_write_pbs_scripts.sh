#!/bin/bash

main() {
  write_cmd "$@"
  arg_parse "$@"
  script=/home/s.benjamin/repo_for_reizel_lab/run_on_zeus/rrbs/bismark_rrbs_pipeline.sh

  if [[ ! -d $OUTPUT_DIR ]]; then
    mkdir $OUTPUT_DIR
  fi
  cd $OUTPUT_DIR

  sample_dirs=$(ls $RAW_SAMPLES_DIR)
  for dir_name in $sample_dirs; do
    if [[ $READ_TYPE == "single_end" ]]; then
      fastq=$(find $RAW_SAMPLES_DIR/$dir_name/ -regex '.*\.\(fastq\|fq\).*')
      script_args=$(echo $extra_trim_galore_opts $extra_meth_extract_opts -genome $genome $non_directional $ovation -n_cores $N_CORES -single-end -input_fastq_file $fastq \> $dir_name.log 2\>\&1)
    else #TODO test paired end, maybe add _1 _2 not only R1 R2
      r1=$(find $RAW_SAMPLES_DIR/$dir_name/ -name R1)
      r2=$(find $RAW_SAMPLES_DIR/$dir_name/ -name R2)
      script_args=$(echo $extra_trim_galore_opts $extra_meth_extract_opts -genome $genome $non_directional $ovation -n_cores $N_CORES -paired-end -paired_input_fastq_files ${r1} ${r2} \> $dir_name.log 2\>\&1)
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

  suggest_submit
}

suggest_submit() {
  echo You may find this list of submission commands usefull:
  find . -name "*.q" | awk '{print "qsub " $1}'
  echo
  printf 'Submit all jobs now? (y/n) '
  read answer
  if [ "$answer" != "${answer#[Yy]}" ]; then # (the #[] operator) means that the variable $answer where any Y or y in 1st position will be dropped if they exist.
    find . -name "*.q" | awk '{print "qsub " $1}' | bash
  fi
  echo Good luck!
}

write_cmd() {
  if [[ $# -gt 2 ]]; then #don't (re)write cmd.txt if no args
    echo \# the command used to prepare the jobs. Note that parentheses are lost >cmd.txt
    echo \# and need to be added to rerun: -extra-trim-galore-options \"multiple quoted options\" >>cmd.txt
    echo "$0" "$@" >>cmd.txt #TODO: preserve quotes that may be in args
  fi
}

arg_parse() {
  [[ $# -eq 0 ]] && {
    help
    exit 1
  }
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
    -ovation)
      ovation="-ovation"
      shift
      ;;
    -extra-trim-galore-options)
      extra_trim_galore_opts=$(echo -extra_trim_galore_opts \'"$2"\')
      shift
      shift
      ;;
    -extra-meth-extractor-options)
      extra_meth_extract_opts=$(echo -extra_meth_extract_opts \'"$2"\')
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
  -genome <mm9/mm10/hg38>

  Optional:
  =========
  [-ovation (must use for Ovation RRBS kit - enables diversity adapters trimming)]
  [-n_cores <int> (DEFAULT: 20)]
  [-non-directional]
  [-extra-meth-extractor-options "multiple quoted arguments"]
  [-extra-trim-galore-options "multiple quoted arguments"]


EOF
}

main "$@"
