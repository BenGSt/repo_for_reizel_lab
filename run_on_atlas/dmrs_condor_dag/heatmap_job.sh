#!/bin/bash 
#usage: <sample_dir>

#first make all_samples_100bp_tiles.bed
	#bedtools unionbedg -names $(find | grep 100bp | grep -v 56pBS_S7_L002_001 | awk -F / '{print $2}' | sed 's/_001//') -header -filler NA -i $(find | grep 100bp | grep -v 
	#56pBS_S7_L002_001) > all_samples_100bp_tiles.bed

source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022
#export PATH=$PATH:/Local/bfe_reizel/anaconda3/envs/ovation_rrbs_pipeline_2022/bin
#TODO: add bedtools to dmrs_pipeline_2022 - 1.7.2022-done?

function main() #args : <path to all_samples_100bp_tiles.bed> <sample_dir>
{
    cd $2
    mkdir heatmaps
    meth_scores_file=$1
    for dmrs in dmrs_25p_hyper.bed dmrs_25p_hypo.bed;  do
        samp_meth_scores=heatmaps/meth_scores_${dmrs}
        head -1 $meth_scores_file > $samp_meth_scores
        bedtools intersect -wa -a $meth_scores_file -b $dmrs >> $samp_meth_scores

	echo make png figures
 	Rscript /srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/make_heatmap.R $3 $4 $5 $6 --scores_bed_file $samp_meth_scores --output_file $(echo $samp_meth_scores| sed 's/.bed/.png/')
    done
}

main `realpath ./all_samples_100bp_tiles.bed`  "$@"
