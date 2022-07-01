#!/usr/bin/bash
source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022
echo Running: Rscript  ~/scripts/repo_for_reizel_lab/scripts_2022/find_dmrs_methylkit.R "$@"
Rscript  ~/scripts/repo_for_reizel_lab/scripts_2022/find_dmrs_methylkit.R "$@"
