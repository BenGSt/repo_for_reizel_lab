#!/bin/bash

N_CORES=10
N_PARALLEL_INSTANCES=4 #each instance uses ~3 cores
BUFFER_SIZE=10G #buffer size for unix sort
MEM=3GB


help()
{
	cat << EOF
	run after trim_illumina_adaptors.sh, trim_diversity_adaptors.sh, align_to_genome.sh
	resources: $N_CORES cores, $MEM RAM

  produces coverage report (only CpG for to include all Cytosines use option --CX)

	-output-dir
	-keep-bam
	-keep-trimmed-fq

	optional:
	-ignore_r2 <int>
	  from Bismark User Guide:
	  ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing results only.
	  Since the first couple of bases in Read 2 of BS-Seq experiments show a severe bias towards non-methylation
	  as a result of end-repairing sonicated fragments with unmethylated cytosines (see M-bias plot),
	  it is recommended that the first couple of bp of Read 2 are removed before starting downstream analysis.
	  Please see the section on M-bias plots in the Bismark User Guide for more details.
EOF
}


main()
{
  arg_parse "$@"
	source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
	cd "$output_dir" || exit 1
	script_name=$(echo $0 | awk -F / '{print $NF}')

	echo
	echo
	echo \################################
	echo \################################
	echo running: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo

	time methylation_calling

	#cleanup
  rm -v $(find ./ | grep -P 'OT|OB')
  rm -v $(find . -name "*.bedGraph.gz")

	echo
	echo
	echo \################################
	echo \################################
	echo finished: $script_name
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo
}


methylation_calling()
{
  alignment_output=$(find . -name '*bismark*deduplicated*bam')
  echo $alignment_output #debug
  echo $alignment_output | grep 'pe' && paired="-p" || paired=""
#  command=$(echo bismark_methylation_extractor --bedgraph $paired $ignore_r2 --multicore $N_PARALLEL_INSTANCES --gzip --buffer_size $BUFFER_SIZE --output methylation_extractor_output $alignment_output)
  command=$(echo bismark_methylation_extractor --bedgraph $paired $ignore_r2 --multicore $N_PARALLEL_INSTANCES --gzip --buffer_size $BUFFER_SIZE $alignment_output)
  echo $SCRIPT_NAME runnig: $command
	$command
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
      -ignore_r2)
        ignore_r2=$(echo --ignore_r2 "$2")
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