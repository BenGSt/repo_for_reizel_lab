#!/usr/bin/env bash

#!/bin/bash

REPO_FOR_REIZEL_LAB=/storage/bfe_reizel/bengst/repo_for_reizel_lab

write_split_job_submission_file() {
  cat <<EOF >condor_submission_files/${sample_name}/split_fastq_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/split_fastq.sh
Arguments = \$(args)
request_cpus = 3
RequestMemory = 250MB
universe = vanilla
log = $(pwd)/logs/$sample_name/${sample_name}_split_fastq.log
output = $(pwd)/logs/$sample_name/${sample_name}_split_fastq.out
error = $(pwd)/logs/$sample_name/${sample_name}_split_fastq.out
queue args from (
  $(
    if [[ $single_end -eq 1 ]]; then
      echo -output-dir $(pwd)/$sample_name/$split/$chunk -chunks $n_chunks -reads-per-chunk $n_reads_per_chunk -input-fastq-file $(realpath $raw_data_dir/$sample_name/*.fastq.gz)
    else
      echo -output-dir $(pwd)/$sample_name/$split/$chunk -chunks $n_chunks -reads-per-chunk $n_reads_per_chunk -paired-input-fastq-files $(realpath $raw_data_dir/$sample_name/*.fastq.gz)
    fi
  )
)
#NOTE: may want to gzip fq files after splitting to save disk space (at the cost of more cpu time)
EOF
}

write_trim_jobs_submission_file() {
  #TODO: if file is not split than no files at $(pwd)/$sample_name/$split/$chunk/\*.fq. make trim jpb use raw_dir or unzip files to $(pwd)/$sample_name/$split/$chunk/\*.fq
  chunk=$1
  if [[ $chunk ]]; then
    filename=condor_submission_files/${sample_name}/trim_job_${sample_name}_${chunk}.sub
  else
    filename=condor_submission_files/${sample_name}/trim_job_${sample_name}.sub
  fi
  cat <<EOF >$filename
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/trim_illumina_adaptors.sh
Arguments = \$(args)
request_cpus = 3
RequestMemory = 500MB
universe = vanilla
log = $(pwd)/logs/$sample_name/\$(name)_trim.log
output = $(pwd)/logs/$sample_name/\$(name)_trim.out
error = $(pwd)/logs/$sample_name/\$(name)_trim.out
queue name, args from (
$(
    if [[ $single_end -eq 1 ]]; then
      echo $sample_name$sep$chunk, \" -output-dir $(pwd)/$sample_name/$split/$chunk -input-fastq-file $(pwd)/$sample_name/$split/$chunk/\*.fq $extra_trim_opts\"
    else
      echo $sample_name$sep$chunk, \" -output-dir $(pwd)/$sample_name/$split/$chunk -paired-input-fastq-files $(pwd)/$sample_name/$split/$chunk/\*.fq $extra_trim_opts\"
    fi
  )
)
#NOTE: If storage turns out to be a bottle neck, may want to gzip fq files after trimming (and / or after splitting)
#      to save disk space (at the cost of more cpu time).

EOF
}

write_align_sub_file() {
  chunk=$1
  if [[ $chunk ]]; then
    filename=condor_submission_files/${sample_name}/bismark_align_job_${sample_name}_${chunk}.sub
  else
    filename=condor_submission_files/${sample_name}/bismark_align_job_${sample_name}.sub
  fi
  cat <<EOF >$filename
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/bismark_align.sh
Arguments = \$(args)
request_cpus = 3
RequestMemory = 40GB
universe = vanilla
log = $(pwd)/logs/$sample_name/\$(name)_bismark_align.log
output = $(pwd)/logs/$sample_name/\$(name)_bismark_align.out
error = $(pwd)/logs/$sample_name/\$(name)_bismark_align.out
queue name, args from (
$(

    if [[ $single_end -eq 1 ]]; then
      echo $sample_name$sep$chunk, -output-dir $(pwd)/$sample_name/$split/$chunk -single-end $non_directional -genome $genome $dovetail
    else
      echo $sample_name$sep$chunk, -output-dir $(pwd)/$sample_name/$split/$chunk -paired-end $non_directional -genome $genome $dovetail
    fi
  )
)
EOF
}

write_deduplicate_job_submission_file() {
  if [[ $single_end -eq 1 ]]; then
    pe_or_se="-single-end"
  else
    pe_or_se="-paired-end"
  fi
  cat <<EOF >condor_submission_files/${sample_name}/deduplicate_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/deduplicate.sh
Arguments = \$(args)
request_cpus = 2
RequestMemory = 25GB
universe = vanilla
log = $(pwd)/logs/$sample_name/\$(name)_deduplicate.log
output = $(pwd)/logs/$sample_name/\$(name)_deduplicate.out
error = $(pwd)/logs/$sample_name/\$(name)_deduplicate.out
queue name, args from (
 $sample_name, $(pwd)/$sample_name $split $pe_or_se
)
EOF
}

write_methylation_calling_job_submission_file() {
  bam_dir=$1 #bam_dir is the directory containing the bam files to be used during m-bias correction (run_mbias.sh)
  cat <<EOF >condor_submission_files/${sample_name}/methylation_calling_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/methylation_calling.sh
Arguments = \$(args)
request_cpus = 3
RequestMemory = 4GB
universe = vanilla
log = $(pwd)/logs/$sample_name/\$(name)_methylation_calling.log
output = $(pwd)/logs/$sample_name/\$(name)_methylation_calling.out
error = $(pwd)/logs/$sample_name/\$(name)_methylation_calling.out
queue name, args from (
  $sample_name, " -output-dir $(pwd)/$sample_name $keep_trimmed_fq $extra_meth_opts $bam_dir"
)
EOF
}

write_bam2nuc_job_submission_file() {
  if [[ $1 == "-override_genome" ]]; then
    genome=$2 #used for mbias correction (fix_mbias.sh)
  fi
  cat <<EOF >condor_submission_files/${sample_name}/bam2nuc_job_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/nucleotide_coverage_report.sh
Arguments = \$(args)
request_cpus = 2
RequestMemory = 10GB
universe = vanilla
log = $(pwd)/logs/$sample_name/\$(name)_bam2nuc.log
output = $(pwd)/logs/$sample_name/\$(name)_bam2nuc.out
error = $(pwd)/logs/$sample_name/\$(name)_bam2nuc.out
queue name, args from (
  $sample_name, -output-dir $(pwd)/$sample_name -genome $genome
)
EOF
}

write_make_tiles_job_submission_file() {
  cat <<EOF >condor_submission_files/${sample_name}/make_tiles_${sample_name}.sub
Initialdir = $(pwd)
executable = $REPO_FOR_REIZEL_LAB/run_on_atlas/bismark_wgbs/make_tiles.sh
Arguments = \$(args)
request_cpus = 1
RequestMemory = 30GB
universe = vanilla
log = $(pwd)/logs/$sample_name/\$(name)_make_tiles.log
output = $(pwd)/logs/$sample_name/\$(name)_make_tiles.out
error = $(pwd)/logs/$sample_name/\$(name)_make_tiles.out
queue name, args from (
  $sample_name, -output-dir $(pwd)/$sample_name -genome $genome
)
EOF
}

write_multiqc_job_submission_file() {
  cat <<EOF >condor_submission_files/multiqc_job.sub
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
}

count_reads() {
  echo "Counting reads in $sample_name to see if the fastq file(s) should be split into chunks"
  n_reads=$(($(pigz -p 1 -cd $(find $raw_data_dir/$sample_name/ -name "*.fastq.gz" | head -1) | wc -l) / 4))
  n_chunks=$((n_reads / n_reads_per_chunk))

  if [[ $((n_reads % n_reads_per_chunk)) -gt 0 ]]; then
    ((n_chunks++)) # add one more chunk for the remainder reads
    n_full_chunks=$((n_chunks - 1))
    remainder_msg=" + 1 chunk of $((n_reads % n_reads_per_chunk)) reads"
  else
    n_full_chunks=$n_chunks
    remainder_msg=
  fi
  echo "n_reads: $n_reads, n_reads_per_chunk: $n_reads_per_chunk"
}

write_split_trim_and_align_sub_files() {
  if [[ $n_reads -gt $n_reads_per_chunk ]]; then
    echo "fastq files will be split into $n_full_chunks chunks of $n_reads_per_chunk reads each" "$remainder_msg"
    echo
    write_split_job_submission_file

    #write condor sub files for jobs to trim and align each chunk
    split="split"
    sep="_"
    for chunk in $(seq -w 00 $((n_chunks - 1))); do
      write_trim_jobs_submission_file $chunk
      write_align_sub_file $chunk
    done
  else # no splitting of fastq files
    write_trim_jobs_submission_file
    write_align_sub_file
  fi
}

write_sample_dag_file() {
  #TODO: refactor this to tidy readable code
  outfile=condor_submission_files/${sample_name}/bismark_wgbs_${sample_name}.dag
  if [[ $split ]]; then
    echo JOB split_fastq $(realpath ./condor_submission_files/${sample_name}/split_fastq_${sample_name}.sub) >$outfile
  else
    truncate --size=0 $outfile #delete previous file's content if it exists
  fi

  if [[ ! $bias_fix ]]; then
    cat <<EOF >>$outfile

$(
      n=0
      for trim_job in $(find ./condor_submission_files/$sample_name/ -name "trim_job_${sample_name}*sub" | sort); do
        echo JOB trim_and_qc_$((n++)) $(realpath $trim_job)
        echo
      done
    )

$(
      n=0
      for align_job in $(find ./condor_submission_files/$sample_name/ -name "bismark_align_job_${sample_name}*sub"| sort); do
        echo JOB bismark_align_$((n++)) $(realpath $align_job)
        echo
      done
      echo $((--n)) >temp_n_value
    )

EOF

    echo JOB deduplicate $(realpath ./condor_submission_files/$sample_name/deduplicate_job_${sample_name}.sub) >>$outfile
  fi
  cat <<EOF >>$outfile

JOB meth_call $(realpath ./condor_submission_files/$sample_name/methylation_calling_job_${sample_name}.sub)

JOB make_tiles $(realpath ./condor_submission_files/$sample_name/make_tiles_${sample_name}.sub)

JOB bam2nuc $(realpath ./condor_submission_files/$sample_name/bam2nuc_job_${sample_name}.sub)


EOF
  if [[ ! $bias_fix ]]; then
    if [[ $split ]]; then
      echo PARENT split_job CHILD \
        $(
          n=$(cat temp_n_value)
          for i in $(seq -w 00 $n); do printf "trim_and_qc_%d " $i; done
        ) >>$outfile
    fi
    cat <<EOF >>$outfile
$(
      n=$(cat temp_n_value)
      rm temp_n_value
      # trim_and_qc -> bismark_align
      printf "PARENT "
      for i in $(seq -w 00 $n); do
        printf "trim_and_qc_%d " $i
      done
      printf "CHILD "
      for i in $(seq -w 00 $n); do
        printf "bismark_align_%d " $i
      done
      printf "\n"

      # bismark_align -> deduplicate
      printf "PARENT "
      for i in $(seq -w 00 $n); do
        printf "bismark_align_%d " $i
      done
      printf "CHILD deduplicate\n"
    )
PARENT deduplicate CHILD meth_call
EOF
  fi
  echo PARENT meth_call CHILD make_tiles bam2nuc >>$outfile
}

write_top_level_dag() {
  mbias_fix=$1
  rm -f ./condor_submission_files/submit_all_bismark_wgbs.dag #incase rerunning the script without delete
  sample_dags=$(realpath $(find ./condor_submission_files/ -name "*.dag" | sort))
  fileout=condor_submission_files/submit_all_bismark_wgbs.dag
  touch $fileout

  if [[ -z $mbias_fix ]]; then
    for sample_name in $(find -L $raw_data_dir -type d | awk -F / 'NR>1{print $NF}' | sort); do
      sample_names+=($sample_name)
    done
  fi

  i=0
  for dag in $sample_dags; do
    echo SUBDAG EXTERNAL ${sample_names[$i]} $dag >>$fileout
    echo PRIORITY ${sample_names[$i]} $i >>$fileout
    echo >>$fileout
    ((i++))
  done
  echo JOB multiqc $(realpath ./condor_submission_files/multiqc_job.sub) >>$fileout
  echo >>$fileout
  #all samples submitted at once
  echo PARENT $(for ((k = 0; k <= $i; k++)); do printf "%s " ${sample_names[$k]}; done) CHILD multiqc >>$fileout
}
