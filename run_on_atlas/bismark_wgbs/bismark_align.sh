#!/bin/bash

source /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/shared.sh



help() {
  cat <<EOF
----------------------------------------
Project: Reizel Lab Bioinformatics Pipelines
Pipeline: Bismark WGBS
Script: bismark_align.sh
Author: Ben G. Steinberg
Last Update: 4 Sep 2023
----------------------------------------

Run after trim_illumina_adaptors.sh

USAGE: bismark_align.sh -single-end or -paired-end -output-dir <path> -genome <mm10 or hg38>
                        [-non-directional] [-keep-trimmed-fq] [-no-dovetail]

Resources: $ALIGN_JOB_CPUS cores, $ALIGN_JOB_MEM RAM

Arguments:
-single-end or -paired-end
-output-dir <path>  Trimmed fq files are expected to be in this directory, output bam will be written here.
-genome <mm10 or hg38>
[-non-directional]   Instructs Bismark to use all four alignment outputs (OT, CTOT, OB, CTOB)
[-keep-trimmed-fq]   Keep the trimmed fq files (default: delete them)
[-no-dovetail]       Don't consider dovetailing alignments as concordant (default: consider them concordant)
                     Explanation from Bismark manual:
                     It is possible, though unusual, for the mates to "dovetail", with the mates seemingly extending
                     "past" each other as in this example:
                         Mate 1:                 GTCAGCTACGATATTGTTTGGGGTGACACATTACGC
                         Mate 2:            TATGAGTCAGCTACGATATTGTTTGGGGTGACACAT
                         Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC
                     Dovetailing is considered inconsistent with concordant alignment, but by default Bismark calls
                     Bowtie 2 with --dovetail, causing it to consider dovetailing alignments as concordant. This
                     becomes relevant whenever reads are clipped from their 5' end prior to mapping, e.g. because of
                     quality or bias issues such as in PBAT or EM-seq libraries.
                     Specify -no-dovetail to turn off this behaviour for paired-end libraries.

EOF
}

main() {
  script_name=$(echo $0 | awk -F / '{print $NF}')
  arg_parse "$@"
  cd "$output_dir" || exit 1
  print_info "running: " "$script_name " "$@"

  echo "cleaning up any previous runs so that restarted jobs aren't corrupted"
  rm -fv $(find . -name "*.bam" -o -name "*unmapped_reads*" -o -name "*ambiguous_reads*" -o -name "*C_to_T*" -o -name "*G_to_A*")
  echo

  #set bismark_genome_location
  if [[ $genome == "mm10" ]]; then
    bismark_genome_location=$MM10_REF # (defined in shared.sh)
  elif [[ $genome == "hg38" ]]; then
    bismark_genome_location=$HG38_REF

  else
    echo genome not recognized
    exit 1
  fi

  align_to_genome

  #cleanup
  rm_fq="rm -v *.fq" #the non gz trimmed fq
  if [[ $keep_trimmed_fq -eq 0 ]]; then
    $rm_fq
  fi
  rm -v *.fq.gz #rm unmapped, ambiguous

  print_info "finished: " "$script_name " "$@"
}

align_to_genome() {
  #see http://felixkrueger.github.io/Bismark/Docs/ :
  #"--parallel 4 for e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~48GB of RAM,
  # but at the same time reduce the alignment time to ~25-30%. You have been warned."
  # Atlas max cpu request is 10 so I want to have 2 instances of bismark (5 cores each theoretically)
  # This is set in align_jobs.sub .

  #fixes Bad file descriptor error (Seems like a bug), and reduces memory usage.
  unmapped_ambig="--un --ambiguous"

  if [[ $read_type == "single_end" ]]; then
    trim_galore_output=$(find . -name '*trimmed.fq*')
    command=$(echo bismark --multicore $BISMARK_INSTANCES --bowtie2 $no_dovetail --genome $bismark_genome_location $trim_galore_output $non_directional $unmapped_ambig)
  else
    trim_galore_output_1=$(find . -name '*val_1.fq*')
    trim_galore_output_2=$(find . -name '*val_2.fq*')
    command=$(echo bismark --multicore $BISMARK_INSTANCES --bowtie2 $no_dovetail --genome $bismark_genome_location -1 $trim_galore_output_1 -2 $trim_galore_output_2 $non_directional $unmapped_ambig)
  fi

  echo runnig: $command
  $command
}

arg_parse() {
  while [[ $# -gt 0 ]]; do
    case $1 in
    -single-end)
      read_type="single_end"
      shift
      ;;
    -paired-end)
      read_type="paired_end"
      shift
      ;;
    -genome)
      genome=$2
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
    -no-dovetail)
      no_dovetail="--no_dovetail"
      shift
      ;;
    -keep-trimmed-fq)
      keep_trimmed_fq=1
      shift
      ;;
    -* | --*)
      help
      exit 1
      ;;
    -h | --help)
      help
      exit 1
      ;;
    esac
  done
}

main "$@"
