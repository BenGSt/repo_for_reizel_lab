#!/bin/bash

main()
{
  if [[ $# -lt 2 ]]; then
    echo USAGE: $0  \<\raw_data_dir\>
    echo raw_data_dir should contain a dir for each sample containing it\'s fastq files.
    echo Run from the directory you wish the output to be written to.
    exit 1
  fi
  write_condor_submition_files $1
  write_condor_dag
  mkdir logs

  echo Submit the jobs by running: condor_submit_dag rna_seq_jobs.dag
  echo Good Luck!
  #TODO: deseq2
  #TODO: single end read option
}

write_condor_submition_files()
{
  raw_dir=$1
  cat << EOF > hisat2_jobs.sub
Initialdir = $(pwd)
executable = hisat2_job.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 8GB
universe = vanilla
log = logs/\$(name)_hisat2.log
output = logs/\$(name)_hisat2.out
error = logs/\$(name)_hisat2.out
queue name, args from (
$( for samp_dir in $(find $raw_dir/* -type d); do echo $samp_dir | awk -F / '{printf $NF", "}'; find $samp_dir | grep -E 'R1|R2'| sort | awk '{printf $0" "}' ; echo $samp_dir | awk -F / '{print "./"$NF"/"$NF".hisat2_output.bam ./"$NF"/"$NF".hisat2.summary.txt "}'  ; done)
)
EOF

  cat << EOF > htseq_jobs.sub
Initialdir = $(pwd)
executable = htseq_job.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 4GB
universe = vanilla
log = logs/\$(name)_htseq-count.log
output = logs/\$(name)_htseq-count.out
error = logs/\$(name)_htseq-count.out
queue name, args from (
$( for samp_dir in $(find $raw_dir/* -type d); do echo $samp_dir | awk -F / '{printf $NF", "}'; echo $samp_dir | awk -F / '{print "./"$NF"/"$NF".hisat2_output.bam ./"$NF"/"$NF".htseq-count_outout.txt"}'  ; done)
)
EOF
}

write_condor_dag()
{
    cat << EOF > rna_seq_jobs.dag
JOB align_hisat2 hisat2_jobs.sub
JOB count_reads_htseq htseq_jobs.sub

PARENT align_hisat2  CHILD count_reads_htseq
EOF
}