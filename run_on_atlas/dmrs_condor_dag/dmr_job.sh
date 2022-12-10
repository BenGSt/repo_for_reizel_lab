#!/usr/bin/bash
source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022
echo Running: Rscript  $REPO_FOR_REIZEL_LAB/scripts_2022/find_dmrs_methylkit.R "$@"
Rscript  $REPO_FOR_REIZEL_LAB/scripts_2022/find_dmrs_methylkit.R "$@"
