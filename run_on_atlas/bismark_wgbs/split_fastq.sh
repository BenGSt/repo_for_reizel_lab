#!/usr/bin/env bash

samp_dir=$1
n_reads_per_file=$(( 4 * $2 )) #4 lines per read
n_chunks=$3
cd $samp_dir || exit 1

# Split fastq files, N_READS_PER_FILE reads per file
mkdir split
for fastq in *.fq; do
  split -dl $n_reads_per_file $fastq split/$(echo $fastq | sed 's/.fastq|fq//')_chunk_ --additional-suffix=.fq
done

cd split
for chunk in  $( seq -w 00 $n_chunks);do
  mkdir $chunk; mv *$chunk.fq $chunk;
done
