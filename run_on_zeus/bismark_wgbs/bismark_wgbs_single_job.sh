#!/bin/bash

N_CORES_TRIM_GALORE=20
N_BISMARK_INSTANCES=4
GENOMIC_REFERENCE_DATA=/home/s.benjamin/genomic_reference_data/

help() {
  cat <<EOF
Author: Ben steinberg 2023
Runs the bismark wgbs pipline on Zeus cluster.
Intended to be used with wgbs_prepare_submission.sh which writes a PBS submission file for each sample, see help there.
EOF
}

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  print_command_info $script_name "$@"
  arg_parse "$@"
  mkdir -p $output_dir
  cd $output_dir || exit 1
  eval "$(micromamba shell hook --shell=bash)"
  micromamba activate /home/s.benjamin/micromamba/envs/wgbs_bismark_pipeline_2023


  if [[ $correct_mbias -eq 1 ]]; then
    methylation_calling
    calculate_tiles 100 10
    write_html_report
  else
    trim_reads_and_fastqc $input_fastq_1 $input_fastq_2
    align_to_genome
    remove_duplicates
    methylation_calling &
    nucleotide_cov_report &
    wait %1 && echo "Done waiting for methylation_calling to complete"
    calculate_tiles 100 10
    (wait %2 && echo "Done waiting for nucleotide_cov_report to complete") || echo "nucleotide_cov_report already done"
    write_html_report
  fi

  echo
  echo
  echo \#################################
  echo \#################################
  echo finished: $script_name "$@"
  echo date: $(date)
  echo \#################################
  echo \#################################
  echo
  echo
}

print_command_info() { # $cmd
  cat <<EOF


#################################
#################################
running: $@
date: $(date)
pwd: $(pwd)
hostname: $(hostname)
#################################
#################################


EOF
}

trim_reads_and_fastqc() { # R1 R2
  # positional argument are  R1, R2 fastq file/s to trim
  # note on multicore from cutadapt manual:
  #   To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
  # note from trim_galore manual:
  #   It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
  #   --cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
  if [[ $read_type == "single_end" ]]; then
    cmd="trim_galore $1 --dont_gzip --cores $N_CORES_TRIM_GALORE --fastqc $extra_trim_galore_opts"
  else #if [[ $read_type == "paired_end" ]]
    cmd="trim_galore --dont_gzip --paired $1 $2 --cores $N_CORES_TRIM_GALORE --fastqc $extra_trim_galore_opts"
  fi

  print_command_info "$cmd"
  $cmd
}

trim_reads() { # R1 R2
  # positional argument are  R1, R2 fastq file/s to trim
  # note on multicore from cutadapt manual:
  #   To automatically detect the number of available cores, use -j 0 (or --cores=0). The detection takes into account resource restrictions that may be in place. For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used. (This works if the cluster systems uses the cpuset(1) mechanism to impose the resource limitation.)
  # note from trim_galore manual:
  #   It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
  #   --cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
  if [[ $read_type == "single_end" ]]; then
    cmd="trim_galore $1 --dont_gzip --cores $N_CORES_TRIM_GALORE $extra_trim_galore_opts"
  else #if [[ $read_type == "paired_end" ]]
    cmd="trim_galore --dont_gzip --paired $1 $2 --cores $N_CORES_TRIM_GALORE  $extra_trim_galore_opts"
  fi

  print_command_info "$cmd"
  $cmd
}

set_bismark_genome_location() {
  if [[ $genome == "mm10" ]]; then
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/from_huji/mm10/Sequence/WholeGenomeFasta/
  elif [[ $genome == "hg38" ]]; then
    #    bismark_genome_location=$GENOMIC_REFERENCE_DATA/hg38/ #Full analysis set
    bismark_genome_location=$GENOMIC_REFERENCE_DATA/hg38/minChromSet/hg38.minChromSet.chroms/
  else
    echo ERROR: genome $genome not supported
    exit 1
  fi
}

align_to_genome() {
  #see http://felixkrueger.github.io/Bismark/Docs/ :
  #"--parallel 4 for e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~48GB of RAM,
  # but at the same time reduce the alignment time to ~25-30%. You have been warned."

  set_bismark_genome_location

  #on Atlas: fixes Bad file descriptor error (Seems like a bug), and reduces memory usage. #TODO: same on Zeus?
  unmapped_ambig="--un --ambiguous"

  if [[ $read_type == "single_end" ]]; then
    trim_galore_output=$(find . -name '*trimmed.fq*')
    command=$(echo bismark --multicore $N_BISMARK_INSTANCES --bowtie2 -genome $bismark_genome_location $trim_galore_output $non_directional $unmapped_ambig)
  else
    trim_galore_output_1=$(find . -name '*val_1.fq*')
    trim_galore_output_2=$(find . -name '*val_2.fq*')
    command=$(echo bismark --multicore $N_BISMARK_INSTANCES --bowtie2 --genome $bismark_genome_location -1 $trim_galore_output_1 -2 $trim_galore_output_2 $non_directional $unmapped_ambig)
  fi

  print_command_info "$command"
  $command || exit 1

  #cleanup
  rm_fq="rm -v *.fq"
  if [[ $keep_trimmed_fq -eq 0 ]]; then
    $rm_fq # rm the non gz trimmed fq
  fi
  rm -v *.fq.gz # rm unmapped, ambiguous
}

remove_duplicates() {
  print_command_info "$(echo deduplicate_bismark ./*bismark*bam)"
  deduplicate_bismark ./*bismark*bam
  rm -v $(find . -name '*.bam' | grep -v deduplicated) #delete bam with duplicates
}

methylation_calling() {
  if [[ $correct_mbias -eq 1 ]]; then
    alignment_output=$(find $biased_dir -name '*bismark*deduplicated*bam')
  else
    alignment_output=$(find . -name '*bismark*deduplicated*bam')
  fi
  echo $alignment_output | grep 'pe' && paired="-p" || paired=""
  # option 1 : use unix sort
  # command=$(echo bismark_methylation_extractor --bedgraph $paired --multicore $N_PARALLEL_INSTANCES --gzip --buffer_size $BUFFER_SIZE $extra_meth_extract_opts $alignment_output)

  #option2 use arrays
  # --ample_memory speeds things up for samples over 10 million reads or so. since it may take over an hour to get going ATLAS policy holds the jobs.
  command=$(echo bismark_methylation_extractor --ample_memory --bedgraph $paired --multicore $N_BISMARK_INSTANCES --gzip $extra_meth_extract_opts $alignment_output)

  print_command_info "$command"
  $command

  #cleanup
  rm -v $(find ./ | grep -P 'OT|OB')
  rm -v $(find . -name "*.bedGraph.gz")
}

nucleotide_cov_report() {
  cmd="$(echo bam2nuc --genome_folder $bismark_genome_location ./*.bam)"
  print_command_info "$cmd"
  $cmd
}

calculate_tiles() {
  # positional args: <tile_size> <min_coverage>
  # adapted from Adam's script.
  # to make WholeGenome_100bpTiles.bed :
  # cat hg38.chrom.sizes  | awk '{OFS="\t"; for (i=1; i<=$2; i+=100) print $1,i,i+100-1}' > hg38_100bp_tiles.bed
  tile_size=$1
  min_coverage=$2

  meth_calling_output=$(find . -name "*.cov.gz")
  print_command_info "calculate_tiles: meth_calling_output=$meth_calling_output"

  if [[ $genome == "mm10" ]]; then
    tiles_file=$GENOMIC_REFERENCE_DATA/from_huji/mm10/mm10_${tile_size}bp_tiles.bed
  elif [[ $genome == "hg38" ]]; then
    tiles_file=$GENOMIC_REFERENCE_DATA/hg38/hg38_${tile_size}bp_tiles.bed
  fi

  # tiles variant 2: First calculate the tiles and then remove tiles with total coverage < 10
  output_file=$(echo ${meth_calling_output} | sed 's/_R[1-2].*//' | awk -v tile_size=$tile_size '{print $1 "_" tile_size "bp_tiles.bed" }')
  bedtools intersect -a $tiles_file -b ${meth_calling_output} -wa -wb | awk -v cov=${min_coverage} -v tileSize=${tile_size} 'BEGIN {OFS="\t"; Prev=-1} {if ($2 == Prev) {T=T+$8+$9; M=M+$8} else {if (Prev!=-1 && T>=cov) {print PrevChr,Prev,Prev+tileSize-1,M/T};T=$8+$9; M=$8;}; Prev=$2; PrevChr=$1}' >${output_file}

  #to unite all tiles from different samples:
  #bedtools unionbedg -names `du -a -L | grep Tiles | awk '{print $2}' | sort | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'` -header -filler NA -i `du -a -L | grep Tiles | awk '{print $2}' | sort` > 100bpTiles_Tiles_Cov10_Tissues.bed
}

write_html_report() {
  cmd=$(echo bismark2report --splitting_report *splitting_report.txt --mbias_report *M-bias.txt \
                            --nucleotide_report *nucleotide_stats.txt --dedup_report *deduplication_report.txt)
  print_command_info "$cmd"
  $cmd

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
    -input-fastq-file)
      read_type="single_end"
      input_fastq_1="$2"
      shift # past argument
      shift # past value
      ;;
    -paired-input-fastq-files)
      read_type="paired_end"
      input_fastq_1="$2"
      shift
      input_fastq_2="$2"
      shift
      shift
      ;;
    -output-dir)
      output_dir="$2"
      shift
      shift
      ;;
    -non-directional)
      non_directional="--non_directional"
      shift
      ;;
    -keep-trimed-fq)
      keep_trimmed_fq=1
      shift
      ;;
    -genome)
      genome=$2
      shift
      shift
      ;;
    -extra-trim-galore-options)
      extra_trim_galore_opts=$2
      shift
      shift
      ;;
    -extra-meth-extractor-options)
      extra_meth_extract_opts=$2
      shift
      shift
      ;;
    -correct-mbias)
      correct_mbias=1
      shift
      ;;
    -biased-dir)
      biased_dir=$2
      shift
      shift
      ;;
    *)
      echo "Unknown option $1"
      help
      exit 1
      ;;
    esac
  done
}

main "$@"
