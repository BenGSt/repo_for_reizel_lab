#!/bin/bash

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab

help()
{
    echo Run The WGBS bismark pipeline:
    echo USAGE: $(echo $0 | awk -F / '{print$NF}') \<1 for single end or 2 for piared end\> \<\raw_data_dir\>
    echo raw_data_dir should contain a dir for each sample containing it\'s fastq files.
    echo Run from the directory you wish the output to be written to.
    echo
    echo Possibly edit the submission files \(you can do this before running the pipeline or after, running additional jobs\).
    echo products: fastqc report, bismark covaregae file, [bam file containing alignments], [100 bp tiles with methylation levels]

}


main()
{


 if [[ $# -lt 2 ]]; then
    help
    exit 1
  fi

  if [[ $1 -eq 1 ]]; then
    single_end=1
  elif [[ $1 -eq 2 ]]; then
    single_end=0
  else
    help
    exit 1
  fi

  write_condor_submition_files $2
  write_condor_dag
  mkdir logs

  echo Submit the jobs by running: condor_submit_dag rna_seq_jobs.dag
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
$(   for sample_name in $(find $raw_dir -type d | awk -F / 'NR>1{print $NF}'|sort) ; do
    if [[ $single_end -eq 1 ]]; then
      echo $sample_name, -output_dir $(pwd)/$sample_name -single-end -input_fastq_file $(realpath $raw_dir/$sample_name/*.fastq.gz)
    else
      echo $sample_name, -output_dir $(pwd)/$sample_name -paired-end -paired_input_fastq_files $(realpath $raw_dir/$sample_name/*.fastq.gz)
    fi
  done
)
)
EOF


  cat << EOF > htseq_jobs.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/htseq_job.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 4GB
universe = vanilla
log = logs/\$(name)_htseq-count.log
output = logs/\$(name)_htseq-count.out
error = logs/\$(name)_htseq-count.out
queue name, args from (
$( for samp_dir in $(find $raw_dir/* -type d); do echo $samp_dir | awk -F / '{printf $NF", "}'; echo $samp_dir | awk -F / '{print "./"$NF"/"$NF".hisat2_output.bam ./"$NF"/"$NF".htseq-count_output.txt"}'  ; done)
)
EOF

#2 conditions for deg, assuming each fastq dir follows the pattern: NAME[0-9].*
condition_1=`ls  $raw_dir/ | sed 's/\(^.*\)[0-9].*$/\1/' | uniq | awk 'NR==1'`
condition_2=`ls  $raw_dir/ | sed 's/\(^.*\)[0-9].*$/\1/' | uniq | awk 'NR==2'`
  cat << EOF > deseq2_job.sub
environment = REPO_FOR_REIZEL_LAB=$REPO_FOR_REIZEL_LAB
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/deseq2_job.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 4GB
universe = vanilla
log = logs/deseq2.log
output = logs/deseq2.out
error = logs/deseq2.out
queue args from (
--htseq_output_dir ./htseq_output --report_dir ./deseq2_deg_results --padj_cutoff 0.01 --log2_fc_cutoff 1 --contrast  c(\"condition\",\"$condition_1\",\"$condition_2\") --csv ./deseq2_results.csv
)
EOF


}}

write_condor_dag()
{
    cat << EOF > rna_seq_jobs.dag
JOB align_hisat2 hisat2_jobs.sub
JOB count_reads_htseq htseq_jobs.sub
JOB find_deg_deseq2 deseq2_job.sub

PARENT align_hisat2  CHILD count_reads_htseq
PARENT count_reads_htseq  CHILD find_deg_deseq2
EOF
}

main "$@"