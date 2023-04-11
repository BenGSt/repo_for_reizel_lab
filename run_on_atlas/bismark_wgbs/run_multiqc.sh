#!/bin/bash

main()
{
  arg_parse "$@"
  if [[ $keep_bam -eq 0 ]]; then
    rm -v $(find . -name "*.bam")
  else
    echo keeping bam files
  fi
  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
  multiqc $multiqc_args
}


arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -h|--help)
        help
        exit 1
        ;;
      -multiqc-args)
        multiqc_args=$2
        shift
        shift
        ;;
      -keep-bam)
        keep_bam=1
        shift
        ;;
      *)
        help
        exit 1
        ;;
    esac
  done
}

help()
{
  cat << EOF
Delete deduplicated bam files and run multiqc on all samples.
-multiqc-args
[-keep-bam]

EOF
}

main "$@"
