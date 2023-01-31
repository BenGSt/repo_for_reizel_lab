#!/bin/bash

help()
{
  cat << EOF
-output-dir <sample_dir>
-genome <hg38 or mm10>

From Bismark Docs:
===================
The script bam2nuc reads BAM files and calculates the mono- and di-nucleotide coverage of the reads
(using the genomic sequence rather than the observed sequence in the reads themselves)
and compares it to the average genomic sequence composition. Reads harbouring InDels are not taken into consideration.
Mono- or dinucleotides containing Ns are ignored as well.
EOF
}


main()
{
	script_name=$(echo $0 | awk -F / '{print $NF}')
	echo
	echo
	echo \#################################
	echo \#################################
	echo running: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo

  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
  arg_parse "$@"
	cd "$output_dir" || exit 1
	mkdir bam2nuc

  if [[ $genome == "mm10" ]]; then
      bismark_genome_location=/storage/bfe_reizel/bengst/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta
    elif [[ $genome == "hg38" ]]; then
      bismark_genome_location=/storage/bfe_reizel/bengst/genomic_reference_data/hg38/analysisSet/hg38.analysisSet.chroms/
    else
      echo genome not recognized
      exit 1
  fi

  bam2nuc --dir bam2nuc --genome_folder $bismark_genome_location ./*.bam
}


arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -h|--help)
        help
        exit 1
        ;;
      -output-dir)
        output_dir="$2"
        shift
        shift
        ;;
      -genome)
        genome=$2
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