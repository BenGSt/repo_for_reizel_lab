#!/bin/bash
# usage: <sample_dir>
source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022

function main()
{
    cd $1
    for dmrs in dmrs_25p_hyper.bed dmrs_25p_hypo.bed; do
        /srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/homer.sh $dmrs $1_dmrs_plus_random_50000_100bp_tiles.bed
    done
}
main "$@"
