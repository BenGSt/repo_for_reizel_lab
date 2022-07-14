#!/bin/bash
# usage: <sample_dir>
source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022

function main()
{
    cd $1
    for dmrs in $(ls| grep -P 'dmrs_[0-9]*p_hyper.bed|dmrs_[0-9]*p_hypo.bed'); do
        /srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/homer.sh $dmrs $1_dmrs_plus_random_50000_100bp_tiles.bed
    done
}
main "$@"
