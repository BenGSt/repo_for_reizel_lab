#!/usr/bin/env bash

#adtaptor trimming
trim_galore -o output --dont_gzip --paired $read1_fq $read2_fq

#alignment
#from https://github.com/ENCODE-DCC/dna-me-pipeline/blob/master/dnanexus/dme-align-pe/resources/usr/bin/dname_align_pe.sh
bismark --bowtie2 -N 1 -L 28 --output_dir output --temp_dir output --multi $ncores \
          input -I $min_insert -X $max_insert -1 ${read1_root}_trimmed.fq -2 ${read2_root}_trimmed.fq
#Standard alignments use a multi-seed length of 20bp with 0 mismatches. These parameters can be modified using the options -L and -N, respectively

#lets try default first
bismark --multicore $N_CORES --genome <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}



#run
/mnt/c/Users/User/Nextcloud/Tzachi_bioinformatics/bioinformatics_repo/run_on_atlas/bismark_wgbs/trim_illumina_adaptors.sh -paired-end -paired_input_fastq_files $(realpath  ../raw_data/100k_reads_PrEC_4/*fastq.gz)    -output_dir 100k_reads_PrEC_4/

/mnt/c/Users/User/Nextcloud/Tzachi_bioinformatics/bioinformatics_repo/run_on_atlas/bismark_wgbs/bismark_align.sh -paired-end  -output-dir ./100k_reads_PrEC_4/
/mnt/c/Users/User/Nextcloud/Tzachi_bioinformatics/bioinformatics_repo/run_on_atlas/bismark_wgbs/deduplicate.sh ./*.bam
/mnt/c/Users/User/Nextcloud/Tzachi_bioinformatics/bioinformatics_repo/run_on_atlas/bismark_wgbs/methylation_calling.sh