#!/bin/bash
#USAGE: $0 <regions_file_mm10.bed> <bg_file.fa>
#results dir in pwd

reference_genome=/utemp/s.benjamin/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta/genome.fa

#to make bg (superset of test regions random tiles):
	#bedtools getfasta -fi $reference_genome -bed "bg_"$len".bed"  -fo "bg_"$len".fa"
	
	#e.g.: bedtools getfasta -fi $reference_genome -bed ./dmrs_plus_random_50000_100bp_tiles.bed  -fo dmrs_plus_random_50000_100bp_tiles.fa^
 
conda activate homer

#Extract Sequence from bed file:
bedtools getfasta -fi $reference_genome -bed $1 -fo `echo $1 | sed s/".bed"/".fa"/`

#find motifs wih rrbs background
 findMotifs.pl `echo $1 | sed s/".bed"/".fa"/` fasta `echo $1 | sed s/".bed"/"_HOMER_RTESULTS"/` -p 20 -fastaBg $2 

conda deactivate