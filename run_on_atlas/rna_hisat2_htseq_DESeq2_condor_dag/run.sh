#!/bin/bash

help() {
  echo Run The RNA-seq deg pipeline:
  echo USAGE: $(echo $0 | awk -F / '{print$NF}') \{-single-end or -paired-end\} -raw-dir \<raw_data_dir\> -genome \<mm10 or hg38\>
  echo raw_data_dir should contain a dir for each sample containing it\'s fastq files.
  echo Biological replicates of the same condition should be in dirs named the same except for tailing numbers e.g. cntrl1 cntrl2.
  echo Run from the directory you wish the output to be written to.
  echo
  echo Possibly edit the submission files \(you can do this before running the pipeline or after, running additional jobs\).
  echo Single end mode isn\'t set up yet. contact me if you need this feature. Ben
  echo products: volcano plot, html report, csv.
  echo
  echo example args pe: -paired-end -raw-dir ~/storage/raw_data/ihep_rna_20.09.2022/X201SC22082747-Z01-F001/01.RawData/ -genome hg38
  echo example args se: -single-end -raw-dir ~/storage/raw_data/rnaseq_row_data_alona/  -genome mm10
  cat << EOF
examples of directory structure:
paired-end with 3 conditions (F, I, ND):
  F1:
  F1_1.fq.gz  F1_2.fq.gz  MD5.txt

  F2:
  F2_1.fq.gz  F2_2.fq.gz  MD5.txt

  F3:
  F3_1.fq.gz  F3_2.fq.gz  MD5.txt

  I4:
  I4_1.fq.gz  I4_2.fq.gz  MD5.txt

  I5:
  I5_1.fq.gz  I5_2.fq.gz  MD5.txt

  I6:
  I6_1.fq.gz  I6_2.fq.gz  MD5.txt

  ND7:
  MD5.txt  ND7_1.fq.gz  ND7_2.fq.gz

  ND8:
  MD5.txt  ND8_1.fq.gz  ND8_2.fq.gz

  ND9:
  MD5.txt  ND9_1.fq.gz  ND9_2.fq.gz

single-end with 2 conditions:
  down1:
  Beta1_down_S11_R1_001.fastq.gz

  down2:
  Beta2_down_S13_R1_001.fastq.gz

  down3:
  Beta3_down_S15_R1_001.fastq.gz

  up1:
  Beta1_up_S10_R1_001.fastq.gz

  up2:
  Beta2_up_S12_R1_001.fastq.gz

  up3:
  Beta3_up_S14_R1_001.fastq.gz
EOF

}

main() {
  REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab
  arg_parse "$@"

  write_condor_submition_files
  write_condor_dag
  mkdir logs

  echo Submit the jobs by running: condor_submit_dag rna_seq_jobs.dag
  echo Good Luck!

}

write_condor_submition_files() {
  cat <<EOF >hisat2_jobs.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/hisat2_job.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 8GB
universe = vanilla
log = logs/\$(name)_hisat2.log
output = logs/\$(name)_hisat2.out
error = logs/\$(name)_hisat2.out
EOF
  if [[ $read_type == "single_end" ]]; then
    cat <<EOF >>hisat2_jobs.sub
    queue name, args from (
$(for samp_dir in $(find $raw_dir/* -type d); do
      echo $samp_dir | awk -F / '{printf $NF", "}'
      printf " -single-end -genome %s" $genome
      printf " -r1 %s " $(find $samp_dir |grep -E "*.fq|*.fq.gz|*.fastq|*.fastq.gz")
      echo $samp_dir | awk -F / '{print "-output-file ./"$NF"/"$NF".hisat2_output.bam -summary-file ./"$NF"/"$NF".hisat2.summary.txt "}'
    done)
)
EOF
  else #if paired end
    cat <<EOF >>hisat2_jobs.sub
    queue name, args from (
$(for samp_dir in $(find $raw_dir/* -type d); do
      r1=$( find $samp_dir | grep -E '.fq.gz|.fastq.gz'| grep -E '_1|R1')
      r2=$( find $samp_dir | grep -E '.fq.gz|.fastq.gz'| grep -E '_2|R2')
      echo $samp_dir | awk -F / '{printf $NF", "}'
      printf " -paired-end -genome %s" $genome
      printf " -r1 %s " $r1
      printf " -r2 %s " $r2
      echo $samp_dir | awk -F / '{print "-output-file ./"$NF"/"$NF".hisat2_output.bam -summary-file ./"$NF"/"$NF".hisat2.summary.txt "}'
    done)
)
EOF
  fi


  cat <<EOF >htseq_jobs.sub
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
$(for samp_dir in $(find $raw_dir/* -type d); do
    echo $samp_dir | awk -F / '{printf $NF", "}'
    echo $samp_dir | awk -v genome=$genome -F / '{print "./"$NF"/"$NF".hisat2_output.bam ./"$NF"/"$NF".htseq-count_output.txt " genome}'
  done)
)
EOF

  #2 conditions for deg, assuming each fastq dir follows the pattern: NAME[0-9].*
  condition_1=$(ls $raw_dir/ | sed 's/\(^.*\)[0-9].*$/\1/' | uniq | awk 'NR==1')
  condition_2=$(ls $raw_dir/ | sed 's/\(^.*\)[0-9].*$/\1/' | uniq | awk 'NR==2')
  cat <<EOF >deseq2_job.sub
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

}

write_condor_dag() {
  cat <<EOF >rna_seq_jobs.dag
JOB align_hisat2 hisat2_jobs.sub
JOB count_reads_htseq htseq_jobs.sub
JOB find_deg_deseq2 deseq2_job.sub

PARENT align_hisat2  CHILD count_reads_htseq
PARENT count_reads_htseq  CHILD find_deg_deseq2
EOF
}

arg_parse() {
  if [[ $# -lt 4 ]]; then
    help
    exit 1
  fi
  while [[ $# -gt 0 ]]; do
    case $1 in
    -h | --help)
      help
      exit 1
      ;;
    -single-end)
      read_type="single_end"
      shift
      ;;
    -paired-end)
      read_type="paired_end"
      shift
      ;;
    -genome)
      genome="$2"
      shift
      shift
      ;;
    -raw-dir)
      raw_dir="$2"
      shift
      shift
      ;;
    *)
      help
      exit 1
      ;;
    esac
  done
}


main "$@"
