#!/usr/bin/env bash

# if enviornment variable DMRS_SCRIPT is not set, set it to the default
DMRS_SCRIPT=${DMRS_SCRIPT:-/mnt/c/Users/User/Nextcloud/Tzachi_bioinformatics/bioinformatics_repo/scripts_2024/find_dmrs_methylkit_v3.R}

help() {
  cat <<EOF
The given arguments will be passed on to the R script that runs the DMRs analysis with
 methylKit in a docker/singularity container. The meth_call_files_dir and output_dir
 will be mounted in the container.

Usage: $0 --meth_call_files_dir <directory>
          --samp_ids <numeric vector>
          --treatments <numeric vector>
          --pipeline <string>
          --output_dir <path>
          --genome <string>
          --meth_difference <number>
          --base_cov <number>
          --tile_cov <number>
          --tile_size <number>
          --filt_hi_perc <number or NULL for no filtering>
          --mc_cores <number>

--meth_call_files_dir  directory where the .cov files are (all will be used)
--samp_ids             double quoted vector with the names of the samples
                           separated by whitespace (must match the order of
                           the .cov files)
--treatments           vector with the condition of each sample
                           (0 or 1) separated by "-" the dmrs are
                           found as the difference between 1 - 0
                           groups (1 - treated , 0 - control)
-p, --pipeline             name of the alignment pipeline, it can be
                           either amp, bismark,bismarkCoverage,
                           bismarkCytosineReport or a list. See
                           methylkit documentation for more details.
--output_dir           directory to save the results in
--genome               mm9, mm10, hg38, etc.
--meth_difference          difference in percent for DMRs, default
                           25% [default: 25]
--base_cov             minimum coverage per CpG, default 1
                           [default: 1]
--tile_cov                 minimum coverage per tile, default 10
                           [default: 10]
--tile_size                tile size, default 100 [default: 100]
--filt_hi_perc         filter out bases with coverage above this
                           percentile (#TODO: not sure
                           this is needed for deduplicated WGBS)
                           [default: 99.9]
--mc_cores             number of cores to use for tileMethylCounts(), unite() and
                           calculateDiffMeth(). must be set to 1 in
                           Windows [default: 1]

Example:
$0 --meth_call_files_dir /data/bismark_trimmed/methylation_coverage \\
   --samp_ids "LNCaP_1 LNCaP_2 LNCaP_3 LNCaP_4 LNCaP_5 PrEC_1 PrEC_2 PrEC_3 PrEC_4" \\
   --treatments "1 1 1 1 1 0 0 0 0" \\
   --pipeline bismarkCoverage \\
   --output_dir /lncap_prec_dmrs \\
   --genome hg19 \\
   --meth_difference 25 \\
   --base_cov 10 \\
   --tile_cov 10 \\
   --tile_size 100 \\
   --filt_hi_perc NULL \\
   --mc.cores 10

Example with directory paths containing spaces:
./find_dmrs.sh \
--meth_call_files_dir "/mnt/c/Users/User/OneDrive - Technion/TEMP_IN_PROGRESS/ARBS_meth_prostate_cancer/data/pidelsy_2016_cell_lines/bismark_trimmed/methylation_coverage" \
  --samp_ids "LNCaP_1 LNCaP_2 LNCaP_3 LNCaP_4 LNCaP_5 PrEC_1 PrEC_2 PrEC_3 PrEC_4" \
    --treatments "1 1 1 1 1 0 0 0 0" \
    --pipeline bismarkCoverage \
    --output_dir "/mnt/c/Users/User/OneDrive - Technion/TEMP_IN_PROGRESS/ARBS_meth_prostate_cancer/analysis/DMRs/LNCaP_vs_PrEC_pidelsy2016" \
    --genome hg19 --meth_difference 25 \
    --base_cov 10 --tile_cov 10  --tile_size 100 \
    --filt_hi_perc NULL --mc.cores 10
EOF
}

set_defaults() {
  # Parse the arguments. set the default values as described in the help.
  # if an argument without a default is missing print a message and then call help()
  meth_call_files_dir=""
  samp_ids=""
  treatments=""
  pipeline=""
  output_dir=""
  genome=""
  meth_difference=25
  base_cov=1
  tile_cov=10
  tile_size=100
  filt_hi_perc=99.9
}

check_missing_args() {
  if [ "$meth_call_files_dir" == "" ]; then
    echo "Missing argument: --meth_call_files_dir"
    help
    exit 1
  fi
  if [ "$samp_ids" == "" ]; then
    echo "Missing argument: --samp_ids"
    help
    exit 1
  fi
  if [ "$treatments" == "" ]; then
    echo "Missing argument: --treatments"
    help
    exit 1
  fi
  if [ "$pipeline" == "" ]; then
    echo "Missing argument: --pipeline"
    help
    exit 1
  fi
  if [ "$output_dir" == "" ]; then
    echo "Missing argument: --output_dir"
    help
    exit 1
  fi
  if [ "$genome" == "" ]; then
    echo "Missing argument: --genome"
    help
    exit 1
  fi

}

arg_parse() {
  while [ "$1" != "" ]; do
    case $1 in
    --meth_call_files_dir)
      shift
      meth_call_files_dir="$1"
      ;;
    --samp_ids)
      shift
      samp_ids="$1"
      ;;
    --treatments)
      shift
      treatments="$1"
      ;;
    --pipeline)
      shift
      pipeline=$1
      ;;
    --output_dir)
      shift
      output_dir="$1"
      ;;
    --genome)
      shift
      genome=$1
      ;;
    --meth_difference)
      shift
      meth_difference=$1
      ;;
    --base_cov)
      shift
      base_cov=$1
      ;;
    --tile_cov)
      shift
      tile_cov=$1
      ;;
    --tile_size)
      shift
      tile_size=$1
      ;;
    --filt_hi_perc)
      shift
      filt_hi_perc=$1
      ;;
    --mc.cores)
      shift
      mc_cores=$1
      ;;
    esac
    shift
  done
}

set_defaults
arg_parse "$@"
check_missing_args

#meth_call_files_dir="/mnt/c/Users/User/OneDrive - Technion/TEMP_IN_PROGRESS/ARBS_meth_prostate_cancer/data/pidelsy_2016_cell_lines/bismark_trimmed/methylation_coverage"
#output_dir="/mnt/c/Users/User/OneDrive - Technion/TEMP_IN_PROGRESS/ARBS_meth_prostate_cancer/analysis/DMRs/LNCaP_vs_PrEC_pidelsy2016"

#tesed and works - run the container and then execute the R script in it
#docker run -p 8787:8787 \
#  -v "$meth_call_files_dir":/meth_call_files_dir \
#  -v "$output_dir":/output_dir \
#  -v "$DMRS_SCRIPT":/find_dmrs.R \
#  methylkit_dmrs:11.8.2024
#
#docker exec reverent_wilbur /find_dmrs.R --meth_call_files_dir /meth_call_files_dir \
#  --samp_ids "LNCaP_1 LNCaP_2 LNCaP_3 LNCaP_4 LNCaP_5 PrEC_1 PrEC_2 PrEC_3 PrEC_4" \
#  --treatments "1 1 1 1 1 0 0 0 0" --pipeline bismarkCoverage \
#  --output_dir /output_dir \
#  --genome hg19 --meth_difference 25 --base_cov 10 --tile_cov 10 \
#  --tile_size 100 --filt_hi_perc NULL --mc.cores 10

#tesed and works - run and remove the container
#docker run --rm \
#  -v "$meth_call_files_dir":/meth_call_files_dir \
#  -v "$output_dir":/output_dir \
#  -v "$DMRS_SCRIPT":/find_dmrs.R \
#  methylkit_dmrs:11.8.2024 /find_dmrs.R --meth_call_files_dir /meth_call_files_dir \
#  --samp_ids "LNCaP_1 LNCaP_2 LNCaP_3 LNCaP_4 LNCaP_5 PrEC_1 PrEC_2 PrEC_3 PrEC_4" \
#  --treatments "1 1 1 1 1 0 0 0 0" --pipeline bismarkCoverage \
#  --output_dir /output_dir \
#  --genome hg19 --meth_difference 25 --base_cov 10 --tile_cov 10 \
#  --tile_size 100 --filt_hi_perc NULL --mc.cores 10

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
singularity exec \
  -B "$meth_call_files_dir":/meth_call_files_dir \
  -B "$output_dir":/output_dir \
  -B "$DMRS_SCRIPT":/find_dmrs.R \
  methylkit_dmrs_11_8_2024.img /find_dmrs.R \
    --meth_call_files_dir /meth_call_files_dir \
    --samp_ids "$samp_ids" \
    --treatments "$treatments" \
    --pipeline "$pipeline" \
    --output_dir /output_dir \
    --genome "$genome" \
    --meth_difference "$meth_difference" \
    --base_cov "$base_cov" \
    --tile_cov "$tile_cov" \
    --tile_size "$tile_size" \
    --filt_hi_perc "$filt_hi_perc" \
    --mc.cores "$mc_cores"

