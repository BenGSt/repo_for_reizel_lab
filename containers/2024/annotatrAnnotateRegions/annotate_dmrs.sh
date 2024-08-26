#!/usr/bin/env bash

# if enviornment variable DMRS_SCRIPT is not set, set it to the default
SCRIPT=${SCRIPT:-/mnt/c/Users/User/Nextcloud/Tzachi_bioinformatics/bioinformatics_repo/scripts_2024/annotate_genomic_regions/annotatr/annotate_dmrs_annotatr.R}

help() {
  cat <<EOF
The given arguments will be passed on to the R script that runs the DMRs  analysis with methylKit in a
docker/singularity container. The given paths will be mounted in the container.

Usage: $0 \\
          --hypo_dmrs <path> \\
          --hyper_dmrs <path> \\
          --all_tiles <path> \\
          --genome <string> \\
          [--internet] \\
          [--print_plots] \\
          [--dont_save_plots] \\
          [--output_dir <path>]

--hypo_dmrs           path to GRanges object saved as rds containing hypomethylated DMRs
--hyper_dmrs          path to GRanges object saved as rds containing hypermethylated DMRs
--all_tiles           path to GRanges object saved as rds containing all tiles (the DMRs are a subset of these)
--genome              genome to use for annotations (available on container: hg19, hg38, mm9, mm10)
--internet            use internet to download more than the basic_genes annotations
--print_plots         print the plots to the console. For use in IDE
--dont_save_plots      don't save the plots as png files
--output_dir          directory to save the plots

Example:
$0 --hypo_dmrs /path/to/hypo_dmrs.rds \\
   --hyper_dmrs /path/to/hyper_dmrs.rds \\
   --all_tiles /path/to/all_tiles.rds \\
   --genome hg19 \\
   --save_plots \\
   --output_dir /path/to/output_dir
EOF
}

check_missing_args() {
  # check if the required arguments are given
  if [ "$hypo_dmrs" == "" ]; then
    echo "Missing argument: --hypo_dmrs"
    help
    exit 1
  fi
  if [ "$hyper_dmrs" == "" ]; then
    echo "Missing argument: --hyper_dmrs"
    help
    exit 1
  fi
  if [ "$all_tiles" == "" ]; then
    echo "Missing argument: --all_tiles"
    help
    exit 1
  fi
  if [ "$genome" == "" ]; then
    echo "Missing argument: --genome"
    help
    exit 1
  fi
  # check if --dont_save_plots is not set and --output_dir is missing
  if [ -z "$dont_save_plots" ] &&
    [ "$output_dir" == "" ]; then
    echo "Missing argument: --output_dir"
    help
    exit 1
  fi
}

arg_parse() {
  while [ "$1" != "" ]; do
    case $1 in
    --hypo_dmrs)
      shift
      hypo_dmrs=$1
      ;;
    --hyper_dmrs)
      shift
      hyper_dmrs=$1
      ;;
    --all_tiles)
      shift
      all_tiles=$1
      ;;
    --genome)
      shift
      genome=$1
      ;;
    --internet)
      internet="--internet"
      ;;
    --print_plots)
      print_plots="--print_plots"
      ;;
    --dont_save_plots)
      dont_save_plots="--dont_save_plots"
      ;;
    --output_dir)
      shift
      output_dir=$1
      ;;
    --help)
      help
      exit
      ;;
    *)
      help
      exit 1
      ;;
    esac
    shift
  done
}

arg_parse "$@"
check_missing_args "$@"

#debug
#docker run -it --rm \
#  -v "$hypo_dmrs":/data/hypo_dmrs \
#  -v "$hyper_dmrs":/data/hyper_dmrs \
#  -v "$all_tiles":/data/all_tiles \
#  -v "$output_dir":/output_dir \
#  -v "$SCRIPT":/annotate_dmrs_annotatr.R \
#  annotatr_annotate_regions:15.8.2024 bash

#- run and remove the container
docker run --rm \
  -v "$hypo_dmrs":/data/hypo_dmrs \
  -v "$hyper_dmrs":/data/hyper_dmrs \
  -v "$all_tiles":/data/all_tiles \
  -v "$output_dir":/output_dir \
  -v "$SCRIPT":/annotate_dmrs_annotatr.R \
  annotatr_annotate_regions:15.8.2024 /annotate_dmrs_annotatr.R\
  --hypo_dmrs /data/hypo_dmrs \
  --hyper_dmrs /data/hyper_dmrs \
  --all_tiles /data/all_tiles \
  --genome "$genome" \
  $internet \
  $print_plots \
  $dont_save_plots \
  --output_dir /output_dir

#run from the sh script args and remove container
#docker run --rm \
#  -v "$meth_call_files_dir":/meth_call_files_dir \
#  -v "$output_dir":/output_dir \
#  -v "$DMRS_SCRIPT":/find_dmrs.R \
#  methylkit_dmrs:11.8.2024 \
#  /find_dmrs.R --meth_call_files_dir /meth_call_files_dir \
#  --samp_ids "$samp_ids" \
#  --treatments "$treatments" \
#  --pipeline "$pipeline" \
#  --output_dir /output_dir \
#  --genome "$genome" \
#  --meth_difference "$meth_difference" \
#  --base_cov "$base_cov" \
#  --tile_cov "$tile_cov" \
#  --tile_size "$tile_size" \
#  --filt_hi_perc "$filt_hi_perc" \
#  --mc.cores "$mc_cores"

# execute the R script in the singularity container
#singularity exec \
#  -B "$meth_call_files_dir":/meth_call_files_dir \
#  -B "$output_dir":/output_dir \
#  -B "$DMRS_SCRIPT":/find_dmrs.R \
#  methylkit_dmrs_11_8_2024.img /find_dmrs.R \
#  --meth_call_files_dir /meth_call_files_dir \
#  --samp_ids "$samp_ids" \
#  --treatments "$treatments" \
#  --pipeline "$pipeline" \
#  --output_dir /output_dir \
#  --genome "$genome" \
#  --meth_difference "$meth_difference" \
#  --base_cov "$base_cov" \
#  --tile_cov "$tile_cov" \
#  --tile_size "$tile_size" \
#  --filt_hi_perc "$filt_hi_perc" \
#  --mc.cores "$mc_cores"
