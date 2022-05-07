#!/usr/bin/env bash

# see notes.txt for more details
main()
{
  arg_parse "$@"
  coefficients_file_1=Thompson_2018_all_CpGs_elastic_net_coefficients_mm10.bed
  coefficients_file_2=Thompson_2018_conserved_CpGs_elastic_net_coefficients_mm10.bed
  coefficients_file_3=Thompson_2018_all_CpGs_ridge_coefficients_mm10.bed
  echo
  echo running hovarth clocks
  echo  method 1 - missing CpGs are given methylation 0
  echo  method 2 - missing CpGs are given methylation scores equal to the sample\'s average meth score
  echo  method 3 - substitute missing CpG meth scores with their 100 bp tile score
  echo

 if [[ $method_ -eq 0 ]]; then
   methods=$(echo 1 2 3)
 else
   methods=$method_
 fi

   for coefficients_file in $coefficients_file_1 $coefficients_file_2 $coefficients_file_3 ; do
       echo $coefficients_file | sed 's/_coefficients_mm10.bed//'

     for method in $methods; do
       echo method=$method
       get_age_predictions_ $coefficients_file $sample_dir
       echo
     done
   done
}


get_age_predictions_() #<coefficients_file>, <sample_dir>
{
   #method 1 - missing CpGs are given methylation 0
   #method 2 - missing CpGs are given methylation scores equal to the sample's average meth score
   #method 3 - substitute missing CpG meth score with there 100 bp tile score ####
   intercept=$(awk 'NR==2{print $2}' $1)

   for sample in $(ls $2); do
     printf $sample"\t"
     if [[ $method -eq 1  ]]; then
       bedtools intersect -wa -wb -f 1 -a $1 -b $2/$sample | awk 'BEGIN{print "#chr \tCpG_coordinate\tcoefficient\t%methylation\tcoeffXmeth"}; {printf("%s\t%d\t%.4f\t%.4f\t%.4f\n", $1,$2,$4,($8 / 100), $4*($8/100))}' | awk -v y=$intercept '{sum+=$5}; END{printf("%.2f\n", y + sum)}'

     elif [[ $method -eq 2  ]]; then
       sample_avg_meth=$( zcat $sample_dir/$sample | awk '{avg+=$4} END{print (avg/NR)/100}')
       bedtools subtract -f 1 -a $coefficients_file -b $sample_dir/$sample | awk -v meth_val=$sample_avg_meth '{OFS="\t"; print $1,$2,$4,meth_val, meth_val*$4}'  > missing_CpGs.bed
       cat <(bedtools intersect -wa -wb -f 1 -a $coefficients_file -b $sample_dir/$sample |  awk 'BEGIN{print "#chr\tCpG_coordinate\tcoefficient\t%methylation\tcoeffXmeth"}; {printf("%s\t%d\t%.4f\t%.4f\t%.4f\n", $1,$2,$4,($8 / 100), $4*($8/100))}' ) missing_CpGs.bed | awk -v y=$intercept '{sum+=$5}; END{printf("%.2f\n", y + sum)}'

     elif [[ $method -eq 3  ]]; then
       barcode=$(echo $sample| sed  's/.*_\([A-Z]\{6\}\)_.*$/\1/')
       sample_column=$(awk -v barcode=$barcode 'NR==1{for (i=0;i<=NF;i++) {if ($i==barcode) print i}}' $all_samp_tiles)
       cut -f 1,2,3,$sample_column $all_samp_tiles  > 100_bp_tiles.bed
       bedtools subtract -f 1 -a $coefficients_file -b $sample_dir/$sample  > missing_CpGs.bed
       bedtools intersect -wa -wb -f 1 -a missing_CpGs.bed -b 100_bp_tiles.bed | awk '{print $1, $2, $4, $8, $4*$8}' > missing_CpGs_scored_by_tile.bed
       cat <(bedtools intersect -wa -wb -f 1 -a $coefficients_file -b $sample_dir/$sample |  awk 'BEGIN{print "#chr\tCpG_coordinate\tcoefficient\t%methylation\tcoeffXmeth"}; {printf("%s\t%d\t%.4f\t%.4f\t%.4f\n", $1,$2,$4,($8 / 100), $4*($8/100))}' ) missing_CpGs_scored_by_tile.bed | awk -v y=$intercept '{sum+=$5}; END{printf("%.2f\n", y + sum)}'
     fi
  done
}


arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -sample_dir)
        sample_dir="$2"
        shift # past argument
        shift # past value
        ;;
      -method)
        method_="$2"
        shift # past argument
        shift # past value
        ;;
      -all_samp_tiles)
        all_samp_tiles="$2"
        shift # past argument
        shift # past value
        ;;
      -h|--help)
        help
        exit 1
        ;;
    esac
  done
}


help()
{
cat << EOF
  By Ben Steinberg (2022).
  runs the hovarth clocks from
    Thompson MJ, Chwiałkowska K, Rubbi L, Lusis AJ, Davis RC, Srivastava A, Korstanje R, Churchill GA, Horvath S, Pellegrini M. A multi-tissue full lifespan epigenetic clock for mice. Aging (Albany NY). 2018; 10:2832-2854. https://doi.org/10.18632/aging.101590
  3 clocks:
    - all_CpGs_elastic_net_coefficients
    - conserved_CpGs_elastic_net_coefficients
    - Thompson_2018_all_CpGs_ridge_coefficients_mm10

  -sample_dir)
      The dir where the bismark .cov files are
      e.g. ../Fah_regeneration/bismkark_meth_extractor_output/
  -method)
      one of 0,1,2,3
      #method 1 - missing CpGs are given methylation 0
      #method 2 - missing CpGs are given methylation scores equal to the sample's average meth score
      #method 3 - substitute missing CpG meth scores with thier 100 bp tile score ####
      0 runs all 3 methods
  -all_samp_tiles)
    bed file with the tiles used for method 3 for all samples
    e.g. ../Fah_regeneration/all_samples_100bp_tiles.bed
  -h|--help)

  example:
    cd /mnt/c/Users/bengs/Nextcloud/Tzachi_bioinformatics/Hovarth_mouse_clocks
     ../bioinformatics_repo/Hovarth_mouse_clocks/run_hovarth_clocks.sh -sample_dir ../Fah_regeneration/bismkark_meth_extractor_output/ -method 0 -all_samp_tiles ../Fah_regeneration/all_samples_100bp_tiles.bed
EOF
}

main "$@"
