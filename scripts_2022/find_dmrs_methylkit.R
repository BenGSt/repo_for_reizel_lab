#!/usr/bin/env Rscript



install_packages = function()
{
  if (!require("argparser", quietly = TRUE))
    install.packages("argparser", dependencies = TRUE)
  
  if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr", dependencies = TRUE)
  
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = TRUE)
  
  BiocManager::install("methylKit", dependencies = TRUE)
  BiocManager::install("GenomicFeatures", dependencies = TRUE)
  install.packages("RMariaDB", dependencies = TRUE) # needed for makeTxDbFromUCSC from "GenomicFeatures"
  BiocManager::install("genomation", dependencies = TRUE)
  BiocManager::install("rGREAT", dependencies = TRUE)
  
}



read_meth_call_files = function(meth_call_files_dir, pipeline_, samp_ids, treatments)
{
  meth_call_files = list.files(path=meth_call_files_dir,
                               pattern="*.cov.gz",
                               full.names=TRUE
  )
  methyl_raw_list = methRead(as.list(meth_call_files),
                             sample.id=as.list(samp_ids),
                             assembly="mm10",
                             pipeline=pipeline_,
                             header=FALSE,
                             treatment=treatments,
                             context="CpG",
                             mincov=10
                             
  )
  
  
  return(methyl_raw_list)
}



filter_bases = function(methyl_raw_list)
{
  # It might be useful to filter samples based on coverage. Particularly, 
  # if our samples are suffering from PCR bias it would be useful to discard bases
  # with very high read coverage. Furthermore, we would also like to discard base
  # that have low read coverage, a high enough read coverage will increase the
  # power of the statistical tests. The code below filters a methylRawList and
  # discards bases that have coverage below 10X and also discards the bases
  # that have more than 99.9th percentile of coverage in each sample.
  
  filtered=filterByCoverage(methyl_raw_list,lo.count=10,lo.perc=NULL,
                            hi.count=NULL,hi.perc=99.9)
  return (filtered)
}


make_tiles = function(meth_call_files_dir, pipeline, samp_ids,
                      treatments)
{
  methyl_raw_list = read_meth_call_files(meth_call_files_dir, pipeline,
                                         samp_ids, treatments)
  
  
  methyl_raw_list = filter_bases(methyl_raw_list)
  # getMethylationStats(methyl_raw_list[[3]],plot=T, both.strands=FALSE)
  # getCoverageStats(methyl_raw_list[[3]],plot=T,both.strands=FALSE)
  
  tiles_raw = tileMethylCounts(methyl_raw_list, win.size=100, step.size=100,
                               mc.cores = 1)
  
  # By default, unite function produces bases/regions covered in all samples.
  # That requirement can be relaxed using ???min.per.group??? option in unite function.
  tiles_raw_Cov10_unite=unite(tiles_raw, destrand=FALSE)
}


#' Finding DMRs with methylKit
#'
#' @param meth_call_files_dir directory where the .cov files are (all will be used)
#' @param samp_ids character vector with the names of the samples (must match the order of the .cov files)
#' @param treatments  vector with the condition of each sample (0 or 1) the dmrs
#'                    are found as the difference between  1 - 0 groups (1 - treated , 0 - control)
#' @param pipeline name of the alignment pipeline, it can be either "amp", "bismark","bismarkCoverage", "bismarkCytosineReport" or a list (default:'amp'). See methylkit documentation for more details.
#' @param output_dir directory to save the results in
#' @param known_genes_file for annotation
#'
#' @return
#' @export methdiff files , bed files, figures
#'
#' @examples   
#' meth_call_files_dir="C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/bismkark_meth_extractor_output"
#' pipeline="bismarkCoverage"
#' samp_ids=c("YoungYoung", "YoungYoung", "OldOld", "YoungYoungProlong", "Young", "YoungYoungProlong", "OldOld", "Old", "Young", "OldOld", "Old")
#' treatments=c(1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0 )
#' setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration")
#' dir.create("figures")
#' setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/figures")
main = function(meth_call_files_dir, samp_ids, treatments, pipeline, output_dir, known_genes_file)
{
  if (! dir.exists(output_dir)){dir.create(output_dir)}
  setwd(output_dir)
  if (! dir.exists("figures")){dir.create("figures")}
  setwd(str_c(output_dir, "/figures"))
  
  tiles_raw_Cov10_unite=make_tiles(meth_call_files_dir, pipeline, samp_ids,
                                   treatments)
  
  png(file="correlation_matrix.png",width=1000,height=1000)
  correlation_matrix = getCorrelation(tiles_raw_Cov10_unite,plot=TRUE)
  dev.off()
  png(file="hierarchical_clustering.png",width=2000,height=1000, res=200)
  hc = clusterSamples(tiles_raw_Cov10_unite, dist="correlation", method="ward",
                      plot=TRUE)
  dev.off()
  png(file="pca.png",width=1500,height=1000,res=100)
  PCASamples(tiles_raw_Cov10_unite)
  dev.off()
  
  
  tiles_raw_Cov10_unite_DMRs = calculateDiffMeth(tiles_raw_Cov10_unite)
  # get hyper methylated bases
  dmrs_25p_hyper = getMethylDiff(tiles_raw_Cov10_unite_DMRs,difference=25,
                                 qvalue=0.01,type="hyper")
  # get hypo methylated bases
  dmrs_25p_hypo = getMethylDiff(tiles_raw_Cov10_unite_DMRs,difference=25,
                                qvalue=0.01,type="hypo")
  # visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
  png("meth_diff_per_chr.png")
  diffMethPerChr(tiles_raw_Cov10_unite_DMRs,plot=TRUE,qvalue.cutoff=0.01,
                 meth.cutoff=25)
  dev.off()
  
  
  #write methDiff files
  write.table(tiles_raw_Cov10_unite_DMRs, str_c(output_dir,"/dmrs.tsv"),sep="\t")
  write.table(dmrs_25p_hyper, str_c(output_dir,"/dmrs_25p_hyper.tsv"),sep="\t")
  write.table(dmrs_25p_hypo, str_c(output_dir,"/dmrs_25p_hypo.tsv"),sep="\t")
  
  #write bed files (only chr start end)
  write.table(getData(dmrs_25p_hyper)[,1:3], str_c(output_dir,"/dmrs_25p_hyper.bed"),sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(getData(dmrs_25p_hypo)[,1:3], str_c(output_dir,"/dmrs_25p_hypo.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(getData(tiles_raw_Cov10_unite)[,1:3], str_c(output_dir,"/all_100bp_tiles_united.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  #bg for graet
  dir_name = str_split(output_dir, "/")[[1]] %>% tail(n=1)
  print(str_c("debug: dir_name = ", dir_name))
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 3000)) %>% unique(), str_c(dir_name,"_dmrs_plus_random_3000_100bp_tiles.bed") ,sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 5000)) %>% unique(), str_c(dir_name, "_dmrs_plus_random_5000_100bp_tiles.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 50000)) %>% unique(), str_c(dir_name,"_dmrs_plus_random_50000_100bp_tiles.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  
  
  if (is.null(known_genes_file))
  {
    #get annotation info
    # TODO: add option to point to a KnownGenes.bed file and don't download every time
    mm10KG_txdb <- makeTxDbFromUCSC(genome="mm10", tablename="knownGene")
    bed_path <- file.path(output_dir, "mm10KnownGenes.bed")
    rtracklayer::export(asBED(mm10KG_txdb), bed_path)
  }
  else 
    bed_path = known_genes_file
  
  
  #annotate DMRs
  gene.obj=readTranscriptFeatures(bed_path)
  dmrs_25p_hyper_annotation=annotateWithGeneParts(as(dmrs_25p_hyper,"GRanges"),gene.obj)
  dmrs_25p_hypo_annotation=annotateWithGeneParts(as(dmrs_25p_hypo,"GRanges"),gene.obj)
  
  
  png(file="hypo_annotationn.png",width=1000,height=1000,res=150)
  plotTargetAnnotation(dmrs_25p_hypo_annotation,precedence=TRUE, main="DMRs 25% hypo annotationn") 
  dev.off()
  
  png(file="hyper_annotationn.png",width=1000,height=1000,res=150)
  plotTargetAnnotation(dmrs_25p_hyper_annotation,precedence=TRUE, main="DMRs 25% hyper annotationn") 
  dev.off()
  
  #Gene Ontology analysis via GREAT
  #TODO: ask tzachi about great params, use bg?
  great_job = submitGreatJob(as(dmrs_25p_hypo, "GRanges"), species = "mm10")
  tb = getEnrichmentTables(great_job, download_by = "tsv")
  
  png(file="RegionGeneAssociationGraphs.png",width=10000,height=2500,res=1000)
  res = plotRegionGeneAssociationGraphs(great_job)
  dev.off()
  
  
}



###__main__##
# install_packages()
# if (!require("argparser", quietly = TRUE))
#   install.packages("argparser")
# 
# if (argv$install-packeges) {install_packages()}

library(argparser, quietly=T)
library(methylKit, quietly=TRUE)
library(GenomicFeatures, quietly=T) # for getting annotation info
library(genomation, quietly=T) #for annotating 
library(rGREAT, quietly=T)
library(dplyr, quietly=T)
library(stringr, quietly=T)

# Create a parser
p <- arg_parser("Find DMRs with methylKit")

# Add command line arguments
p <- add_argument(p, "--meth_call_files_dir",  help="directory where the .cov files are (all will be used)", short="-m")
p <- add_argument(p, "--samp_ids", help="vector with the names of the samples (must match the order of the .cov files)", short="-s")
p <- add_argument(p, "--treatments", help="vector with the condition of each sample (0 or 1) the dmrs are found as the difference between  1 - 0 groups (1 - treated , 0 - control)", short="-t")
p <- add_argument(p, "--pipeline", help="name of the alignment pipeline, it can be either amp, bismark,bismarkCoverage, bismarkCytosineReport or a list. See methylkit documentation for more details.", short="-p")
p <- add_argument(p, "--output_dir", help="directory to save the results in", short="-o")
p <- add_argument(p, "--known_genes_file", help="annotaion info e.g. mm10KnownGenes.bed, if none is given will be downloaded")
p <- add_argument(p, "--install-packeges", help="install requirements")

# Parse the command line arguments
argv <- parse_args(p)


treatments = strsplit(argv$treatments,' +')[[1]] %>% as.numeric
samp_ids = strsplit(argv$samp_ids,' +')[[1]]

main(argv$meth_call_files_dir, samp_ids, treatments, argv$pipeline, argv$output_dir, argv$known_genes_file)

