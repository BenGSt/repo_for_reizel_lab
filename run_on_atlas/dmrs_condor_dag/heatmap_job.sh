#!/bin/bash 
#usage: <sample_dir>

#first make all_samples_100bp_tiles.bed
#	bedtools unionbedg -names $(find | grep 100bp | grep -v dmrs | awk -F / '{print $3"_"$4}' | sed 's/\(FGC18.6\)_.*/\1/' | awk -F _ '{OFS="_"; print $1,$3}' | sed 's/_FGC1866//' | sed 's/FGC1876/2/') -header -filler NA -i $(find | grep 100bp | grep -v dmrs) > ./dmrs_01.07.2022/all_samples_100bp_tiles_each_run_as_separate_samples.bed
#TODO: make 100 bp tiles automatically

source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022


function main() #args : <path to all_samples_100bp_tiles.bed> <sample_dir - where output of dmr_job is>
{
    cd $2
    mkdir heatmaps
    meth_scores_file=$1
    for dmrs in $(ls| grep -P 'dmrs_[0-9]*p_hyper.bed|dmrs_[0-9]*p_hypo.bed');  do
        samp_meth_scores=heatmaps/meth_scores_${dmrs}
        head -1 $meth_scores_file > $samp_meth_scores
        bedtools intersect -wa -a $meth_scores_file -b $dmrs >> $samp_meth_scores
        #TODO: this intersection ends up with less tiles than in $dmrs.
        # is it because of the discarding of tiles with less than 10X cov?
        # it shouldn't be because find_dmrs_methylkit.R is supposed to filter all CpGs with less than 10x.
        # what 100bp tiles show up in $dmrs but not in $meth_scores_file ?

	      echo make png figures
 	      Rscript /srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/make_heatmap.R $3 $4 $5 $6 --scores_bed_file $samp_meth_scores --output_file $(echo $samp_meth_scores| sed 's/.bed/.png/')
    done
}

main   "$@"
#TODO: add arg_parse arg path_to_all_samples_100bp_tiles
