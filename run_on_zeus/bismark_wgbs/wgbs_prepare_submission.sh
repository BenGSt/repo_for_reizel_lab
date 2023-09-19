#!/bin/bash

N_CORES=20
REPO_FOR_REIZEL_LAB=/home/s.benjamin/repo_for_reizel_lab

help() {
  cat <<EOF
Usage:  $(echo "$0" | awk -F / '{print$NF}') {-single-end or -paired-end} -raw-data-dir <raw_data_dir>
    -genome <mm10 or hg38>  [optional]

Run from the directory you wish the output to be written to.
raw_data_dir should contain a dir for each sample containing it's fastq files.

Note about methylation bias correction: I recommend running the pipeline once without additional options, you should
then view the m-bias plots in the MultiQC report. The expected unbiased result is a uniform distribution of the
average methylation levels across read positions. If the results are biased, fix this by either running the methylation
calling jobs again ignoring the biased bases (run this script in a different directory with -correct-mbias & -biased-dir
& -extra-meth-extract-options), or running the pipeline again with trimmed reads. Each of these approaches has it's
advantages and disadvantages. Ignoring aligned bases is faster. Trimming the reads may improve alignment, but can use a
lot of resources and time.

optional:
-non-directional
  Use for non directional libraries. Instructs Bismark to align to OT, CTOT, OB, CTOB.

-keep-bam
  Don't delete the deduplicated bam files. Useful for running methylation calling jobs again to fix m-bias without
  trimming and rerunning the pipeline, and possibly other downstream analysis.

-correct-mbias & -biased-dir & -extra-meth-extract-options
  Use these three options together to correct m-bias by rerunning the methylation calling step and ignoring biased
  positions.


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

EOF

}
main() {
  if [[ $# -gt 2 ]]; then #don't (re)write cmd.txt if no args
    echo \# the command used to prepare the jobs. Note that parentheses are lost >cmd.txt
    echo \# and need to be added to rerun: -extra-trim-galore-options \"multiple quoted options\" >>cmd.txt
    echo "$0" "$@" >>cmd.txt #TODO: preserve quotes that may be in args
  fi

  arg_parse "$@"
  script=$REPO_FOR_REIZEL_LAB/run_on_zeus/bismark_wgbs/bismark_wgbs_single_job.sh

  for sample_name in $(find -L $raw_data_dir -type d | awk -F / 'NR>1{print $NF}' | sort); do
    input_fastq=$(realpath $raw_data_dir/$sample_name/*.fastq.gz | tr '\n' ' ')
    if [[ $single_end -eq 1 ]]; then
      args="$correct_mbias $biased_dir -output-dir $(realpath $PWD)/$sample_name -input-fastq-file $input_fastq -genome $genome $non_directional $extra_trim_opts $extra_meth_opts"
    else
      args="$correct_mbias $biased_dir -output-dir $(realpath $PWD)/$sample_name -paired-input-fastq-files $input_fastq -genome $genome $non_directional $extra_trim_opts $extra_meth_opts"
    fi

    mkdir -p $sample_name
    cd $sample_name
    #TODO: 8.5.23: use zeus_new_q, holding all jobs for some reason. using long_q for now
    #TODO: 18.9.2023 : switched back to zeus_new_q. testing
    cat <<EOF >bismark_wgbs_${sample_name}.q
#!/bin/bash
#PBS  -N  bismark_wgbs_${sample_name}
#PBS  -q  zeus_new_q
#PBS  -l select=1:ncpus=${N_CORES}
#PBS  -l select=mem=64gb
PBS_O_WORKDIR=$(realpath $PWD)
cd \$PBS_O_WORKDIR


$script $args > bismark_wgbs_${sample_name}.log 2>&1

EOF
    cd ..
  done

  printf 'Submit all jobs now? (y/n) '
  read answer
  if [ "$answer" != "${answer#[Yy]}" ]; then # this grammar (the #[] operator) means that the variable $answer where any Y or y in 1st position will be dropped if they exist.
    find . -name "*.q" | awk '{print "qsub " $1}' | bash
  else
    echo You may find this list of submission commands usefull:
    find . -name "*.q" | awk '{print "qsub " $1}'
  fi
  echo Good luck!

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
      non_directional="-non_directional"
      shift
      ;;
    -raw-data-dir)
      raw_data_dir=$2
      shift
      shift
      ;;
      #    -keep-bam)
      #      keep_bam="-keep-bam" #TODO: option to keep or delete bam
      #      shift
      #      ;;
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
    -extra-meth-extract-options)
      extra_meth_opts=$(echo -extra-meth_extract-options \'"$2"\')
      shift
      shift
      ;;
    -correct-mbias)
      correct_mbias="-correct-mbias"
      shift
      ;;
    -biased-dir)
      biased_dir="-biased-dir $2"
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
