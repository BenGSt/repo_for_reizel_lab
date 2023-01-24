#!/bin/bash

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab

help()
{
    echo Run The WGBS bismark pipeline:
    echo USAGE: $(echo $0 | awk -F / '{print$NF}') \{-single-end or -paired-end\} -raw-data-dir \<\raw_data_dir\> -genome \<mm10 or hg38\>\[-non-directional\]
    echo
    echo raw_data_dir should contain a dir for each sample containing it\'s fastq files.
    echo -non-directional instructs Bismark to use all four alignment outputs \(OT, CTOT, OB, CTOB\).
    echo Run from the directory you wish the output to be written to.
    echo
    echo Possibly edit the submission files \(you can do this before running the pipeline or after, running additional jobs\).
    echo products: fastqc report, bismark covaregae file, [bam file containing alignments], [100 bp tiles with methylation levels]

}


main()
{
  arg_parse "$@"
  write_condor_submition_files $raw_data_dir
  write_condor_dag
  mkdir logs

  echo Submit the jobs by running: condor_submit_dag bismark_wgbs.dag
  echo Good Luck!
  #TODO: single end read option

}

write_condor_submition_files() # <raw_dir>
{
  raw_dir=$1
  cat << EOF > trim_jobs.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/trim_illumina_adaptors.sh
Arguments = \$(args)
request_cpus = 8
RequestMemory = 500MB
universe = vanilla
log = logs/\$(name)_trim.log
output = logs/\$(name)_trim.out
error = logs/\$(name)_trim.out
queue name, args from (
$(
  for sample_name in $(find $raw_dir -type d | awk -F / 'NR>1{print $NF}'|sort) ; do
    if [[ $single_end -eq 1 ]]; then
      echo $sample_name, -output-dir $(pwd)/$sample_name -single-end -input_fastq_file $(realpath $raw_dir/$sample_name/*.fastq.gz)
    else
      echo $sample_name, -output-dir $(pwd)/$sample_name -paired-end -paired_input_fastq_files $(realpath $raw_dir/$sample_name/*.fastq.gz)
    fi
  done
)
)
EOF


  cat << EOF > bismark_align_jobs.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/bismark_align.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 30GB
universe = vanilla
log = logs/\$(name)_bismark_align.log
output = logs/\$(name)_bismark_align.out
error = logs/\$(name)_bismark_align.out
queue name, args from (
$(
  for sample_name in $(find $raw_dir -type d | awk -F / 'NR>1{print $NF}'|sort) ; do
      if [[ $single_end -eq 1 ]]; then
        echo $sample_name, -output-dir $(pwd)/$sample_name -single-end $non_directional -genome $genome
      else
        echo $sample_name, -output-dir $(pwd)/$sample_name -paired-end $non_directional -genome $genome
      fi
  done
)
)
EOF


  cat << EOF > deduplicate_jobs.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/deduplicate.sh
Arguments = \$(args)
request_cpus = 2
RequestMemory = 20GB
universe = vanilla
log = logs/\$(name)_deduplicate.log
output = logs/\$(name)_deduplicate.out
error = logs/\$(name)_deduplicate.out
queue name, args from (
$(
  for sample_name in $(find $raw_dir -type d | awk -F / 'NR>1{print $NF}'|sort) ; do
        echo $sample_name, $(pwd)/$sample_name
  done
)
)
EOF

  cat << EOF > methylation_calling_jobs.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/methylation_calling.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 3GB
universe = vanilla
log = logs/\$(name)_methylation_calling.log
output = logs/\$(name)_methylation_calling.out
error = logs/\$(name)_methylation_calling.out
queue name, args from (
$(
  for sample_name in $(find $raw_dir -type d | awk -F / 'NR>1{print $NF}'|sort) ; do
        echo $sample_name, -output-dir $(pwd)/$sample_name $keep_bam $keep_trimmed_fq
  done
)
)
EOF

  cat << EOF > make_tiles.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/make_tiles.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 3GB
universe = vanilla
log = logs/\$(name)_make_tiles.log
output = logs/\$(name)_make_tiles.out
error = logs/\$(name)_make_tiles.out
queue name, args from (
$(
  for sample_name in $(find $raw_dir -type d | awk -F / 'NR>1{print $NF}'|sort) ; do
        echo $sample_name, -output-dir $(pwd)/$sample_name -genome $genome
  done
)
)
EOF


}

write_condor_dag()
{
    cat << EOF > bismark_wgbs.dag
JOB trim_and_qc trim_jobs.sub
JOB bismark_align bismark_align_jobs.sub
JOB deduplicate deduplicate_jobs.sub
JOB meth_call methylation_calling_jobs.sub
JOB make_tiles make_tiles.sub

PARENT trim_and_qc  CHILD bismark_align
PARENT bismark_align  CHILD deduplicate
PARENT deduplicate  CHILD meth_call
PARENT meth_call  CHILD make_tiles
EOF
}

arg_parse()
{
  if [[ $# -eq 0 ]]; then
    help
    exit 1
  fi
  while [[ $# -gt 0 ]]; do
    case $1 in
     -single-end)
        single_end=1
        shift
        ;;
     -paired-end)
        single_end=0
        shift
        ;;
 	  -non-directional)
        non_directional="--non_directional"
        shift
        ;;
 	  -raw-data-dir)
        raw_data_dir=$2
        shift
        shift
        ;;
    -keep-bam)
        keep_bam="-keep-bam"
        shift
        ;;
    -keep-trimmed-fq)
        keep_trimmed_fq="-keep-trimmed-fq"
        shift
        ;;
    -genome)
        genome=$2
        shift
        shift
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