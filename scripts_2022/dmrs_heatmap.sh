#!/bin/bash

#first make all_samples_100bp_tiles.bed
	#bedtools unionbedg -names $(find | grep 100bp | grep -v 56pBS_S7_L002_001 | awk -F / '{print $2}' | sed 's/_001//')  -header -filler NA -i $(find | grep 100bp | grep -v 56pBS_S7_L002_001)  > all_samples_100bp_tiles.bed


function main() #args : <path to all_samples_100bp_tiles.bed> 
{
    mkdir heatmaps
    meth_scores_file=$1
    for dmrs in dmrs_25p_hyper.bed dmrs_25p_hypo.bed; do
        outfile=heatmaps/meth_scores_${dmrs}
        head -1 $meth_scores_file > $outfile
        bedtools intersect -wa -a $meth_scores_file -b $dmrs >> $outfile
		
		#make png figures
		Rscript make_heatmap.R --scores_bed_file $outfile  $(echo $outfile| sed 's/.bed/.png/')
    done
}

main "$@"