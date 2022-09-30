#!/bin/bash
#USAGE: $0 <regions_file_mm10.bed> <bg_file.bed>
#results dir in pwd
N_CORES=10
reference_genome=/storage/bfe_reizel/bengst/genomic_reference_data/from_huji/mm10/Sequence/WholeGenomeFasta/genome.fa

#make bg (superset of test regions random tiles):
bedtools getfasta -fi $reference_genome -bed $2  -fo $(echo $2 | sed 's/.bed/.fa/')
	#e.g.: bedtools getfasta -fi $reference_genome -bed ./dmrs_plus_random_50000_100bp_tiles.bed  -fo dmrs_plus_random_50000_100bp_tiles.fa^

#Extract Sequence from bed file:
bedtools getfasta -fi $reference_genome -bed $1 -fo `echo $1 | sed s/".bed"/".fa"/`

#find motifs with  background
findMotifs.pl `echo $1 | sed s/".bed"/".fa"/` fasta `echo $1 | sed s/".bed"/"_HOMER_RESULTS"/` -p $N_CORES -fastaBg $(echo $2 | sed 's/.bed/.fa/')
