#!/bin/bash

help()
{
    echo Run The RNA-seq deg pipeline:
    echo USAGE: $(echo $0 | awk -F / '{print$NF}') \<1 for single end or 2 for piared end\> \<\raw_data_dir\>
    echo raw_data_dir should contain a dir for each sample containing it\'s fastq files.
    echo Biological replicates of the same condition should be in dirs named the same except for tailing numbers e.g. cntrl1 cntrl2.
    echo Run from the directory you wish the output to be written to.
    echo
    echo Possibly edit the submission files \(you can do this before running the pipeline or after, running additional jobs\).
    echo Single end mode isn\'t set up yet. contact me if you need this feature. Ben
    echo products: volcano plot, html report, csv.

}


main()
{
 if [[ $# -lt 2 ]]; then
    help
    exit 1
  fi

  if [[ $1 -eq 1 ]]; then
    single_end=1
    echo single end #debug
  elif [[ $1 -eq 2 ]]; then
    single_end=0
    echo paired_end #debug
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

write_condor_submition_files()
{
  raw_dir=$1
  cat << EOF > hisat2_jobs.sub
Initialdir = $(pwd)
executable = /srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/hisat2_job.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 8GB
universe = vanilla
log = logs/\$(name)_hisat2.log
output = logs/\$(name)_hisat2.out
error = logs/\$(name)_hisat2.out
EOF
#  if [[ $single_end -eq 0 ]]; then
#    cat << EOF >> hisat2_jobs.sub
#    queue name, args from (
#$( for samp_dir in $(find $raw_dir/* -type d); do echo $samp_dir | awk -F / '{printf $NF", "}'; find $samp_dir | grep -E '_1|_2|R1|R2'| sort | awk '{printf $0" "}' ; echo $samp_dir | awk -F / '{print "./"$NF"/"$NF".hisat2_output.bam ./"$NF"/"$NF".hisat2.summary.txt "}'  ; done)
#)
#EOF
#  else #if single end
    cat << EOF >> hisat2_jobs.sub
    queue name, args from (
$( for samp_dir in $(find $raw_dir/* -type d); do echo $samp_dir | awk -F / '{printf $NF", "}'; find $samp_dir | grep -E '.fq.gz|.fastq.gz'| sort | awk '{printf $0" "}' ; echo $samp_dir | awk -F / '{print "./"$NF"/"$NF".hisat2_output.bam ./"$NF"/"$NF".hisat2.summary.txt "}'  ; done)
)
EOF
#  fi

  cat << EOF > htseq_jobs.sub
Initialdir = $(pwd)
executable = /srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/htseq_job.sh
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
Initialdir = $(pwd)
executable = /srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/deseq2_job.sh
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

}

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