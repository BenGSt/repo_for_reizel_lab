#!/usr/bin/env bash

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh
n_reads_per_chunk=25000000 #default value (may be overwritten by arg_parse)

help() {
  cat <<EOF
----------------------------------------
Project: Reizel Lab Bioinformatics Pipelines
Pipeline: Bismark WGBS
Script: prepare_jobs.sh
Author: Ben G. Steinberg
Last Update: 4 Sep 2023
----------------------------------------

This script is the starting point for running the bismark WGBS pipeline on the Atlas cluster. It is designed to be run 3
times for each sample batch, in other words it has 3 modes of operation:
(1) First run, called directly, will write the condor submission file prep.sub, and prompt the user to submit it. A file
    named prep2.cmd be written, which should be used to run the script again after the jobs in prep.sub are done.
    (run "bash prep2.cmd")
(2) Second run, submitted for each sample by prep.sub with the -job flag, the job will count the number of reads in the
    fastq files to determine whether to split them. It will then write the condor submission files
    needed for the sample.
(3) Third run, called with the -top-level option when the user runs "bash prep2.cmd". It will write the submission file
    for MultiQC as well as the top level dag needed to run the batch and prompt the user to submit it.


USAGE: prepare_jobs.sh {-single-end or -paired-end} -raw-data-dir <raw_data_dir> -genome <mm10 or hg38> [optional]

  <raw_data_dir> should contain a dir for each sample containing it's fastq files.
  Run from the directory you wish the output to be written to.


optional arguments:
-n-reads-per-chunk
  Number of reads per chunk for splitting fastq files. Default is $n_reads_per_chunk

-non-directional
  Use for non directional libraries. Instructs Bismark to align to OT, CTOT, OB, CTOB.

-delete-bam
  Delete the deduplicated bam files. Default is to keep them for running methylation calling jobs again to fix m-bias without
  trimming and rerunning the pipeline, and possibly other downstream analysis. If not running methylation calling jobs again,
  bam files should be deleted because they large and not needed for most downstream analysis (use the .cov files).

-extra-meth-extract-options "multiple quoted options"
handy options (from Bismark manual):
=====================================

Ignore bases in aligned reads.
------------------------------------------------------------------------------------------------------------------
--ignore <int>
    Ignore the first <int> bp from the 5' end of Read 1 (or single-end alignment files) when processing
    the methylation call string. This can remove e.g. a restriction enzyme site at the start of each read or any other
    source of bias (such as PBAT-Seq data).

--ignore_r2 <int>
    Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing results only. Since the first couple of
    bases in Read 2 of BS-Seq experiments show a severe bias towards non-methylation as a result of end-repairing
    sonicated fragments with unmethylated cytosines (see M-bias plot), it is recommended that the first couple of
    bp of Read 2 are removed before starting downstream analysis. Please see the section on M-bias plots in the Bismark
    User Guide for more details.

--ignore_3prime <int>
    Ignore the last <int> bp from the 3' end of Read 1 (or single-end alignment files) when processing the methylation
    call string. This can remove unwanted biases from the end of reads.

--ignore_3prime_r2 <int>
    Ignore the last <int> bp from the 3' end of Read 2 of paired-end sequencing results only. This can remove unwanted
    biases from the end of reads.

Other
------------------------------------------------------------------------------------------------------------------------
--no_overlap
    For paired-end reads it is theoretically possible that Read 1 and Read 2 overlap. This option avoids scoring
    overlapping methylation calls twice (only methylation calls of read 1 are used for in the process since read 1 has
    historically higher quality basecalls than read 2). Whilst this option removes a bias towards more methylation calls
    in the center of sequenced fragments it may de facto remove a sizeable proportion of the data. This option is on by
    default for paired-end data but can be disabled using --include_overlap. Default: ON.

--include_overlap
    For paired-end data all methylation calls will be extracted irrespective of whether they overlap or not.
    Default: OFF.

--zero_based
    Write out an additional coverage file (ending in .zero.cov) that uses 0-based genomic start and 1-based genomic end
    coordinates (zero-based, half-open), like used in the bedGraph file, instead of using 1-based coordinates
    throughout. Default: OFF.


-extra-trim-galore-options "multiple quoted options"
handy options (from trim_galore manual):
=====================================

Remove bases from reads before alignment.
------------------------------------------------------------------------------------------------------------------
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


A note about methylation bias correction: After running the pipeline once view the m-bias plots in the MultiQC report.
The expected unbiased result is a uniform distribution of the average methylation levels across read positions. If the
results are biased, fix this by one of 2 methods: (a) use the provided correct_mbias.sh script to the pipeline again
starting from the methylation calling and ignoring the biased bases of the aligned reads. Or (b) running the pipeline
again, trimming the biased bases and re aligning. Each of these approaches has it's advantages and disadvantages.
Ignoring aligned bases is faster, and should be stable when running on Atlas. Trimming the reads may improve alignment,
hence give more accurate results. However, the memory usage may go above 40GB leading to jobs failing
(tested before splitting fastq files), This has yet to be tested with split fastq files. You may want to consider
trimming R1 and R2 symmetrically and/or using the "--dovetail" bismark option for the bowtie2 aligner
(--dovetail is actually the default).

EOF
}

prepare_sample() {
  mkdir -p condor_submission_files/$sample_name
  mkdir -p logs/$sample_name

  # count reads to see if fastq file is longer than n_reads_per_chunk reads, if so it will be split into chunks of
  # length n_reads_per_chunk, given as an argument to this script, defaults to 25M.
  count_reads

  #write sub files for trimming and aligning each chunk (or one sub file if no splitting)
  write_split_trim_and_align_sub_files

  # the following sub files are not dependent on splitting
  write_deduplicate_job_submission_file
  write_methylation_calling_job_submission_file
  write_bam2nuc_job_submission_file
  write_make_tiles_job_submission_file
  write_sample_dag_file
}

build_args_str() {
  if [[ $single_end -eq 1 ]]; then
    args_for_perp_sub="-single-end"
  else
     args_for_perp_sub="-paired-end"
  fi

  args_for_perp_sub="$args_for_perp_sub
    $non_directional \
    -raw-data-dir $raw_data_dir \
    $keep_bam \
    $keep_trimmed_fq \
    $dovetail \
    -genome $genome \
    -n-reads-per-chunk $n_reads_per_chunk \
    $extra_trim_opts \
    ${extra_meth_opts/-extra-options/-extra-meth-extract-options}"
    echo $args_for_perp_sub
}


write_prep_submission_files() {
  samples=$(find -L $raw_data_dir -type d | awk -F / 'NR>1{print $NF}' | sort)
  if [[ -z $samples ]]; then
    echo "ERROR: no samples found in $raw_data_dir"
    exit 1
  fi


  args="\"$(build_args_str) -sample-name \$(sample_name) -job\""
  cat <<EOF >condor_submission_files/prep/prep.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/prepare_jobs.sh
Arguments = $args
request_cpus = 1
RequestMemory = 500MB
universe = vanilla
log = $(pwd)/logs/prep/\$(sample_name).log
output = $(pwd)/logs/prep/\$(sample_name).out
error = $(pwd)/logs/prep/\$(sample_name).out
queue sample_name from(
$(
    for sample_name in $samples; do
      echo $sample_name
    done
  )
)
EOF
}

submit_prep_jobs() {
  echo The following samples will be prepared:
  for sample_name in $samples; do
    echo $sample_name
  done
  echo
  echo To submit the jobs, run: condor_submit condor_submission_files/prep/prep.sub
  echo !NOTE: After the initial prep jobs are finished, run \"bash prep2.cmd\" to prepare and submit the top level dag jobs!
  echo
  printf 'Submit prep jobs now? (y/n) '
  read answer
  if [ "$answer" != "${answer#[Yy]}" ]; then # this grammar (the #[] operator) means that the variable $answer where any Y or y in 1st position will be dropped if they exist.
    condor_submit condor_submission_files/prep/prep.sub
  fi
}

submit_top_level_dag() {
  #TODO: tee output to reminder.txt because instruction are complicated and need to be carried out over long time periods (>= days
  cat <<EOF | tee reminder.txt
To submit the samples sepratly, you may run the following commands:
$(
for dag in $(find ./condor_submission_files/ -name "*.dag" | grep -v submit_all | sort -V); do
    echo condor_submit_dag $dag
done
  )

To submit the top level dag (all samples), run the following commands:
condor_submit_dag condor_submission_files/submit_all_bismark_wgbs.dag

!Note: after the top level dag jobs are finished, correct m-bias by running correct_mbias.sh!
!e.g. correct_mbias.sh --biased_dir \$PWD --ignore 8 --ignore_r2 9 --ignore_3prime_r2 1!
!for more info see $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/correct_mbias.sh --help!

!Note: Atlas policy is to hold jobs after 3 days of running, so you may need to release them using condor_release!

I have seen a case where memory usage in alignment jobs was less than 20% of the requested memory after 1 hour,
and jobs were held. I'm guessing this had to do with io bottleneck in nfs and is a low probability situation.
to deal with this you can wait for the dag to finish (with errors) and then run the rescue dag:
condor_submit_dag -DoRescueFrom 1 condor_submission_files/submit_all_bismark_wgbs.dag

Check the status of the jobs by running: condor_q -nobatch -dag
!Good luck and happy clustering!

EOF
  printf 'Submit top level dag now? (y/n) '
  read answer
  if [ "$answer" != "${answer#[Yy]}" ]; then # this grammar (the #[] operator) means that the variable $answer where any Y or y in 1st position will be dropped if they exist.
    condor_submit_dag condor_submission_files/submit_all_bismark_wgbs.dag
  fi

  #TODO: add instructions for restarting held jobs (draft below), print this whole block with cat <<EOF
  #TODO: update 3.9.2023: may not need this, released jobs restart and delete previous output.
  #TODO: Also, still remains to be seen if the cmd below work.
#  cat <<EOF
#One way to restart a job that is part of a DAG is to use the condor_dagman tool. You can use the -f option
#with condor_dagman to force it to re-run a specific node in the DAG. For example, if your DAG is defined in a
#file called mydag.dag and the job you want to restart is represented by the node NODE_A, you can use the following
#command to force condor_dagman to re-run that node:
#condor_dagman -f NODE_A mydag.dag
#EOF
}

arg_parse() {
  if [[ $# -lt 2 ]]; then
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
    -dovetail) #seems this is on by default.
      dovetail="-dovetail"
      shift
      ;;
    -genome)
      genome=$2
      shift
      shift
      ;;
    -n-reads-per-chunk) #for splitting fastq files, default is 100M
      n_reads_per_chunk=$2
      shift
      shift
      ;;
    -extra-trim-galore-options)
      extra_trim_opts=$(echo -extra-trim-galore-options \'"$2"\')
      shift
      shift
      ;;
    -extra-meth-extract-options)
      extra_meth_opts=$(echo -extra-options \'"$2"\')
      shift
      shift
      ;;
    -job)
      job=1
      shift
      ;;
    -top-level)
      top_level=1
      shift
      ;;
    -sample-name)
      sample_name=$2
      shift
      shift
      ;;
    *)
      echo "ERROR: unknown argument $1"
      help
      exit 1
      ;;
    esac
  done
}

main() {
  mkdir -p logs/prep
  mkdir -p condor_submission_files/prep/
  arg_parse "$@"

  if [[ $job ]]; then
    prepare_sample
  elif [[ $top_level ]]; then
    write_multiqc_job_submission_file
    write_top_level_dag
    submit_top_level_dag
  else
    echo "# -extar-trim-galore-options and -extra-meth-extract-options must be quoted strings, these are lost here. add to run again" >cmd.txt
    echo "$0" "$@" >>cmd.txt
    echo "$0 $(build_args_str) -top-level" >prep2.cmd
    write_prep_submission_files "$@"
    submit_prep_jobs
  fi

}

main "$@"
