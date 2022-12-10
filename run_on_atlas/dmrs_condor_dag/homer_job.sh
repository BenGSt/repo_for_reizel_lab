#!/bin/bash
# usage: <sample_dir>
source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022

function main()
{
    cd $1
    # make random 50k tiles from all_100bp_tiles.bed (better for homer as opposed to great which requires superset)
    shuf -n 50000 all_100bp_tiles_united.bed > random_50000_100bp_tiles_for_homer_bg.bed
    for dmrs in $(ls| grep -P 'dmrs_[0-9]*p_hyper.bed|dmrs_[0-9]*p_hypo.bed'); do
        # old version using great bg - remove after testing - 30.09.2022
        #/srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/homer.sh $dmrs $1_dmrs_plus_random_50000_100bp_tiles.bed
       $REPO_FOR_REIZEL_LAB/scripts_2022/homer.sh $dmrs random_50000_100bp_tiles_for_homer_bg.bed
    done
}
main "$@"
