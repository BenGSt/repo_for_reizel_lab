#!/bin/bash 
#usage: <sample_dir>

#first make all_samples_100bp_tiles.bed
	#bedtools unionbedg -names $(find | grep 100bp | grep -v dmrs | awk -F / '{print $3"_"$4}' | sed 's/\(FGC18.6\)_.*/\1/' | awk -F _ '{OFS="_"; print $1,$3}' | sed 's/_FGC1866//' | sed 's/FGC1876/2/') -header -filler NA -i $(find | grep 100bp | grep -v dmrs) > all_samples_100bp_tiles_each_run_as_separate_samples.bed
#TODO: make 100 bp tiles automatically

source /Local/bfe_reizel/anaconda3/bin/activate dmrs_pipeline_2022


function main() #args : <path to all_samples_100bp_tiles.bed> <sample_dir>
{
    cd $2
    mkdir heatmaps
    meth_scores_file=$1
    for dmrs in dmrs_25p_hyper.bed dmrs_25p_hypo.bed;  do
        samp_meth_scores=heatmaps/meth_scores_${dmrs}
        head -1 $meth_scores_file > $samp_meth_scores
        bedtools intersect -wa -a $meth_scores_file -b $dmrs >> $samp_meth_scores
        #TODO: I'm making 100bp tile scroes twice. Once in the rrbs pipline and once again in find_dmrs_methylkit.R
        #    choose one for better efficiency

	echo make png figures
 	Rscript /srv01/technion/bengst/scripts/repo_for_reizel_lab/scripts_2022/make_heatmap.R $3 $4 $5 $6 --scores_bed_file $samp_meth_scores --output_file $(echo $samp_meth_scores| sed 's/.bed/.png/')
    done
}

main   "$@"
#TODO: add arg_parse arg path_to_all_samples_100bp_tiles
