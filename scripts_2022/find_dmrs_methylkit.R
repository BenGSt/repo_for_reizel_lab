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
                               pattern="*.cov.gz|*.cov",
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

  return(tiles_raw_Cov10_unite)
}


write_meth_scores = function(methylBase.obj, output_file )
{
    table = percMethylation(methylBase.obj, rowids=TRUE)
    rownames(table) = str_replace_all(rownames(table), '\\.', '\t')
    write.table(paste("#chr\tstart\tend", paste(colnames(table),collapse="\t"),sep="\t"),
                output_file, sep="\t", row.names=FALSE , col.names=FALSE,
                quote=FALSE)
    write.table(table / 100, output_file, append=TRUE, sep="\t",
                row.names=TRUE, col.names=FALSE, quote=FALSE)
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
#' @param meth_difference difference in percent for DMRs, default 25%
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
main = function(meth_call_files_dir, samp_ids, treatments, pipeline, output_dir, known_genes_file, meth_difference)
{
  tiles_raw_Cov10_unite=make_tiles(meth_call_files_dir, pipeline, samp_ids,
                                   treatments)
  
  if (! dir.exists(output_dir))
    dir.create(output_dir)
  setwd(output_dir)
  if (! dir.exists("figures"))
    dir.create("figures")
  setwd("./figures")
  
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

  if (is.null(meth_difference))
    meth_difference=25
  # get hyper methylated bases
  dmrs_hyper = getMethylDiff(tiles_raw_Cov10_unite_DMRs,difference=meth_difference,
                                 qvalue=0.01,type="hyper")
  # get hypo methylated bases
  dmrs_hypo = getMethylDiff(tiles_raw_Cov10_unite_DMRs,difference=meth_difference,
                                qvalue=0.01,type="hypo")
  # visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
  png("meth_diff_per_chr.png")
  diffMethPerChr(tiles_raw_Cov10_unite_DMRs,plot=TRUE,qvalue.cutoff=0.01,
                 meth.cutoff=meth_difference)
  dev.off()

  #write methDiff files
  write.table(tiles_raw_Cov10_unite_DMRs, str_c(output_dir,"/dmrs.tsv"),sep="\t")
  write.table(dmrs_hyper, str_c(output_dir, "/dmrs_", meth_difference, "p_hyper.tsv"),sep="\t")
  write.table(dmrs_hypo, str_c(output_dir, "/dmrs_", meth_difference, "p_hypo.tsv"),sep="\t")
  
  #write bed files (only chr start end)
  write.table(getData(dmrs_hyper)[,1:3], str_c(output_dir, "/dmrs_", meth_difference, "p_hyper.bed"),sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(getData(dmrs_hypo)[,1:3], str_c(output_dir, "/dmrs_", meth_difference, "p_hypo.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(getData(tiles_raw_Cov10_unite)[,1:3], str_c(output_dir,"/all_100bp_tiles_united.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)

  #write all 100bp tiles meth scores (not only dmrs) that can be used for heatmaps and other applications instead of
  # the 100bp tiles produced by the rrbs pipeline
  write_meth_scores(tiles_raw_Cov10_unite, str_c(output_dir,"/all_samps_100bp_tiles_meth_scores.bed"))

  #bg for great
  name = str_split(output_dir, "/")[[1]] %>% tail(n=1)
  write.table(rbind(getData(dmrs_hyper)[,1:3], getData(dmrs_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 3000)) %>% unique(), str_c(output_dir,"/", name,"_dmrs_plus_random_3000_100bp_tiles.bed") ,sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_hyper)[,1:3], getData(dmrs_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 5000)) %>% unique(), str_c(output_dir,"/", name, "_dmrs_plus_random_5000_100bp_tiles.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_hyper)[,1:3], getData(dmrs_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 50000)) %>% unique(), str_c(output_dir,"/", name,"_dmrs_plus_random_50000_100bp_tiles.bed"),sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)

  if (is.null(known_genes_file))
  {
    #get annotation info
    #download  KnownGenes.bed file it wasn't given by user
    mm10KG_txdb <- makeTxDbFromUCSC(genome="mm10", tablename="knownGene")
    bed_path <- file.path(output_dir, "mm10KnownGenes.bed")
    rtracklayer::export(asBED(mm10KG_txdb), bed_path)
  }
  else 
    bed_path = known_genes_file

  #annotate DMRs
  gene.obj=readTranscriptFeatures(bed_path)
  dmrs_hyper_annotation=annotateWithGeneParts(as(dmrs_hyper,"GRanges"),gene.obj)
  dmrs_hypo_annotation=annotateWithGeneParts(as(dmrs_hypo,"GRanges"),gene.obj)

  png(file="hypo_annotation.png",width=1000,height=1000,res=150)
  plotTargetAnnotation(dmrs_hypo_annotation,precedence=TRUE, main=str_c("DMRs ", meth_difference, "% hypo annotation"))
  dev.off()
  
  png(file="hyper_annotation.png",width=1000,height=1000,res=150)
  plotTargetAnnotation(dmrs_hyper_annotation,precedence=TRUE, main=str_c("DMRs ", meth_difference, "% hyper annotation"))
  dev.off()
  
  #Gene Ontology analysis via GREAT
  #TODO: ask tzachi about great params, use bg?
  great_job = submitGreatJob(as(dmrs_hypo, "GRanges"), species = "mm10")
  tb = getEnrichmentTables(great_job, download_by = "tsv")
  
  png(file="RegionGeneAssociationGraphs.png",width=10000,height=2500,res=1000)
  res = plotRegionGeneAssociationGraphs(great_job)
  dev.off()
}


suppressMessages(library(argparser))
suppressMessages(library(methylKit))
suppressMessages(library(GenomicFeatures)) # for getting annotation info
suppressMessages(library(genomation)) #for annotating
suppressMessages(library(rGREAT))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# Create a parser
p <- arg_parser("Find DMRs with methylKit")
p <- add_argument(p, "--meth_call_files_dir",  help="directory where the .cov files are (all will be used)", short="-m")
p <- add_argument(p, "--samp_ids", help="vector with the names of the samples separated by \"-\" (must match the order of the .cov files)", short="-s")
p <- add_argument(p, "--treatments", help="vector with the condition of each sample (0 or 1) separated by \"-\" the dmrs are found as the difference between  1 - 0 groups (1 - treated , 0 - control)", short="-t")
p <- add_argument(p, "--pipeline", help="name of the alignment pipeline, it can be either amp, bismark,bismarkCoverage, bismarkCytosineReport or a list. See methylkit documentation for more details.", short="-p")
p <- add_argument(p, "--output_dir", help="directory to save the results in", short="-o")
p <- add_argument(p, "--known_genes_file", help="annotation info e.g. mm10KnownGenes.bed, if none is given will be downloaded")
p <- add_argument(p, "--meth_difference", help="difference in percent for DMRs, default 25%")
#TODO: p <- add_argument(p, "--install-packages", help="install requirements")
argv <- parse_args(p)

# to solve problem on condor multiple jobs will be using "-" instead of whitespace
# to separate treatments and samp_ids
treatments = strsplit(argv$treatments,'-')[[1]] %>% as.numeric
samp_ids = strsplit(argv$samp_ids,'-')[[1]]
# treatments = strsplit(argv$treatments,' +')[[1]] %>% as.numeric
# samp_ids = strsplit(argv$samp_ids,' +')[[1]]

#allow list(...) as pipline input
if (str_detect(argv$pipeline,"list"))
    argv$pipeline = eval(parse(text=argv$pipeline))

main(normalizePath(argv$meth_call_files_dir), samp_ids, treatments, argv$pipeline, normalizePath(argv$output_dir), argv$known_genes_file,
     as.numeric(argv$meth_difference))
#TODO: test if normalizePath() works well