#!/bin/bash

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab

help() {
  echo Run The WGBS bismark pipeline \(separate dag for each sample\):
  echo USAGE: "$(echo "$0" | awk -F / '{print$NF}')" \{-single-end or -paired-end\} -raw-data-dir \<raw_data_dir\> \
    -genome \<mm10 or hg38\>\[-non-directional\] \[-ignore_r2 \<int\>\] \[-extra-trim-galore-options "multiple quoted options"\]
  echo
  echo raw_data_dir should contain a dir for each sample containing it\'s fastq files.
  echo -non-directional instructs Bismark to use all four alignment outputs \(OT, CTOT, OB, CTOB\).
  echo Run from the directory you wish the output to be written to.
  echo
  echo Possibly edit the submission files \(you can do this before running the pipeline or after, running additional jobs\).
  echo products: fastqc report, bismark covaregae file, [bam file containing alignments], [100 bp tiles with methylation levels]
  cat <<EOF
optional:

-ignore_r2 <int>
  from Bismark User Guide:
  ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing results only.
  Since the first couple of bases in Read 2 of BS-Seq experiments show a severe bias towards non-methylation
  as a result of end-repairing sonicated fragments with unmethylated cytosines (see M-bias plot),
  it is recommended that the first couple of bp of Read 2 are removed before starting downstream analysis.
  Please see the section on M-bias plots in the Bismark User Guide for more details.

-extra-trim-galore-options "multiple quoted options"
handy extra options from trim_galore manual:
=====================================
--clip_R1 <int>         Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end
                      reads). This may be useful if the qualities were very poor, or if there is some
                      sort of unwanted bias at the 5' end. Default: OFF.

--clip_R2 <int>         Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end reads
                        only). This may be useful if the qualities were very poor, or if there is some sort
                        of unwanted bias at the 5' end. For paired-end BS-Seq, it is recommended to remove
                        the first few bp because the end-repair reaction may introduce a bias towards low
                        methylation. Please refer to the M-bias plot section in the Bismark User Guide for
                        some examples. Default: OFF.

--three_prime_clip_R1 <int>     Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end
                        reads) AFTER adapter/quality trimming has been performed. This may remove some unwanted
                        bias from the 3' end that is not directly related to adapter sequence or basecall quality.
                        Default: OFF.

--three_prime_clip_R2 <int>     Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER
                        adapter/quality trimming has been performed. This may remove some unwanted bias from
                        the 3' end that is not directly related to adapter sequence or basecall quality.
                        Default: OFF.

EOF

}

main() {
  echo \# the command used to prepare the jobs. Note that parentheses are lost > cmd.txt
  echo \# and need to be added to rerun: -extra-trim-galore-options \"multiple quoted options\" >> cmd.txt
  echo "$0" "$@" >> cmd.txt #TODO: preserve quotes that may be in args
  arg_parse "$@"
  write_condor_submission_files $raw_data_dir
  mkdir logs

  echo Submit the jobs by running: condor_submit_dag ./condor_submission_files/submit_all_bismark_wgbs.dag
  echo Good Luck!
  #TODO: single end read option

}

write_condor_submission_files() { # <raw_dir>
  raw_dir=$1
  sample_names=()
  mkdir condor_submission_files
  for sample_name in $(find -L $raw_dir -type d | awk -F / 'NR>1{print $NF}' | sort); do
    sample_names+=($sample_name)
    cat << EOF >condor_submission_files/trim_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/trim_illumina_adaptors.sh
Arguments = \$(args)
request_cpus = 8
RequestMemory = 500MB
universe = vanilla
log = $(pwd)/logs/\$(name)_trim.log
output = $(pwd)/logs/\$(name)_trim.out
error = $(pwd)/logs/\$(name)_trim.out
queue name, args from (
$(
      if [[ $single_end -eq 1 ]]; then
        echo $sample_name, \" -output-dir $(pwd)/$sample_name -input-fastq-file $(realpath $raw_dir/$sample_name/*.fastq.gz) $extra_trim_opts\"
      else
        echo $sample_name, \" -output-dir $(pwd)/$sample_name -paired-input-fastq-files $(realpath $raw_dir/$sample_name/*.fastq.gz) $extra_trim_opts\"
      fi
    )
)
EOF

    cat << EOF >condor_submission_files/bismark_align_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/bismark_align.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 40GB
universe = vanilla
log = $(pwd)/logs/\$(name)_bismark_align.log
output = $(pwd)/logs/\$(name)_bismark_align.out
error = $(pwd)/logs/\$(name)_bismark_align.out
queue name, args from (
$(

      if [[ $single_end -eq 1 ]]; then
        echo $sample_name, -output-dir $(pwd)/$sample_name -single-end $non_directional -genome $genome
      else
        echo $sample_name, -output-dir $(pwd)/$sample_name -paired-end $non_directional -genome $genome
      fi
    )
)
EOF

    cat << EOF >condor_submission_files/deduplicate_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/deduplicate.sh
Arguments = \$(args)
request_cpus = 2
RequestMemory = 20GB
universe = vanilla
log = $(pwd)/logs/\$(name)_deduplicate.log
output = $(pwd)/logs/\$(name)_deduplicate.out
error = $(pwd)/logs/\$(name)_deduplicate.out
queue name, args from (
 $sample_name, $(pwd)/$sample_name
)
EOF

    cat << EOF >condor_submission_files/methylation_calling_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/methylation_calling.sh
Arguments = \$(args)
request_cpus = 10
RequestMemory = 3GB
universe = vanilla
log = $(pwd)/logs/\$(name)_methylation_calling.log
output = $(pwd)/logs/\$(name)_methylation_calling.out
error = $(pwd)/logs/\$(name)_methylation_calling.out
queue name, args from (
  $sample_name, -output-dir $(pwd)/$sample_name $ignore_r2 $keep_trimmed_fq

)
EOF

    cat << EOF >condor_submission_files/bam2nuc_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/nucleotide_coverage_report.sh
Arguments = \$(args)
request_cpus = 2
RequestMemory = 10GB
universe = vanilla
log = $(pwd)/logs/\$(name)_bam2nuc.log
output = $(pwd)/logs/\$(name)_bam2nuc.out
error = $(pwd)/logs/\$(name)_bam2nuc.out
queue name, args from (
  $sample_name, -output-dir $(pwd)/$sample_name -genome $genome
)
EOF

    cat << EOF >condor_submission_files/make_tiles_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/make_tiles.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 30GB
universe = vanilla
log = $(pwd)/logs/\$(name)_make_tiles.log
output = $(pwd)/logs/\$(name)_make_tiles.out
error = $(pwd)/logs/\$(name)_make_tiles.out
queue name, args from (
  $sample_name, -output-dir $(pwd)/$sample_name -genome $genome
)
EOF


  cat << EOF >condor_submission_files/bismark_wgbs_${sample_name}.dag
JOB trim_and_qc $(realpath ./condor_submission_files/trim_job_${sample_name}.sub)
JOB bismark_align $(realpath ./condor_submission_files/bismark_align_job_${sample_name}.sub)
JOB deduplicate $(realpath ./condor_submission_files/deduplicate_job_${sample_name}.sub)
JOB meth_call $(realpath ./condor_submission_files/methylation_calling_job_${sample_name}.sub)
JOB make_tiles $(realpath ./condor_submission_files/make_tiles_${sample_name}.sub)
JOB bam2nuc $(realpath ./condor_submission_files/bam2nuc_job_${sample_name}.sub)

PARENT trim_and_qc  CHILD bismark_align
PARENT bismark_align  CHILD deduplicate
PARENT deduplicate  CHILD meth_call bam2nuc
PARENT meth_call  CHILD make_tiles
EOF

  done

  cat << EOF > condor_submission_files/multiqc_job.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/run_multiqc.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 500MB
universe = vanilla
log = $(pwd)/logs/multiqc_job.log
output = $(pwd)/logs/multiqc_job.out
error = $(pwd)/logs/multiqc_job.out
queue args from (
  "$keep_bam -multiqc-args '$(pwd) --outdir multiqc'"
)
EOF

  rm -f ./condor_submission_files/submit_all_bismark_wgbs.dag #incase rerunning the script without delete
  sample_dags=$(realpath condor_submission_files/*.dag)
  touch ./condor_submission_files/submit_all_bismark_wgbs.dag
  i=0
  for dag in $sample_dags; do
    echo SUBDAG EXTERNAL ${sample_names[$i]} $dag >> condor_submission_files/submit_all_bismark_wgbs.dag
    echo PRIORITY ${sample_names[$i]} $i >> condor_submission_files/submit_all_bismark_wgbs.dag #TODO: as of 20.3.23 still needs testing.
    echo >> condor_submission_files/submit_all_bismark_wgbs.dag
    ((i++))
  done
  echo JOB multiqc $(realpath ./condor_submission_files/multiqc_job.sub) >> condor_submission_files/submit_all_bismark_wgbs.dag
  echo PARENT $(for ((k=0; k<=$i; k++)); do printf "%s " ${sample_names[$k]}; done) CHILD multiqc >> condor_submission_files/submit_all_bismark_wgbs.dag

}


arg_parse() {
  if [[ $# -eq 0 ]]; then
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
    -extra-trim-galore-options)
      extra_trim_opts=$(echo -extra-trim-galore-options \'"$2"\')
      shift
      shift
      ;;
    -ignore_r2)
      ignore_r2=$(echo -ignore_r2 "$2")
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
