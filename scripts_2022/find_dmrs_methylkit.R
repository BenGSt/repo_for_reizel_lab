library(methylKit)
library(GenomicFeatures) # for getting annotation info
library(genomation) #for annotating 
library(rGREAT)
library(dplyr)

install_packages = function()
{
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("methylKit")
  BiocManager::install("GenomicFeatures")
  install.packages("RMariaDB") # needed for makeTxDbFromUCSC from "GenomicFeatures"
  BiocManager::install("genomation")
  BiocManager::install("rGREAT")
  
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
  # That requirement can be relaxed using “min.per.group” option in unite function.
  tiles_raw_Cov10_unite=unite(tiles_raw, destrand=FALSE)
}


main = function()
{
  meth_call_files_dir="C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/bismkark_meth_extractor_output"
  meth_call_files_dir_no_prolong="C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/bismkark_meth_extractor_output_no_prolong"
  pipeline="bismarkCoverage"
  samp_ids_no_prolong=c("YoungYoung1", "YoungYoung2", "OldOld1", 
             "Young1", "OldOld2", "Old1", "Young2",
             "OldOld3", "Old2")
  treatments_no_prolong=c(1, 1, 1, 0, 1, 0, 0, 1, 0 )
  
  samp_ids=c("YoungYoung", "YoungYoung", "OldOld", "YoungYoungProlong", "Young", "YoungYoungProlong", "OldOld", "Old", "Young", "OldOld", "Old")
  treatments=c(1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0 )
  

  
  setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration")
  dir.create("figures")
  setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/figures")
  
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
  #
  # get hypo methylated bases
  dmrs_25p_hypo = getMethylDiff(tiles_raw_Cov10_unite_DMRs,difference=25,
                                qvalue=0.01,type="hypo")
  
  # visualize the distribution of hypo/hyper-methylated bases/regions per chromosome 
  diffMethPerChr(tiles_raw_Cov10_unite_DMRs,plot=TRUE,qvalue.cutoff=0.01,
                 meth.cutoff=25)
  
  #write methDiff files
  write.table(tiles_raw_Cov10_unite_DMRs, "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs.tsv",sep="\t")
  write.table(dmrs_25p_hyper, "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_25p_hyper.tsv",sep="\t")
  write.table(dmrs_25p_hypo, "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_25p_hypo.tsv",sep="\t")
  
  #write bed files (only chr start end)
  write.table(getData(dmrs_25p_hyper)[,1:3], "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_25p_hyper.bed",sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(getData(dmrs_25p_hypo)[,1:3], "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_25p_hypo.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(getData(tiles_raw_Cov10_unite)[,1:3], "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/all_100bp_tiles_united.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  #bg for graet
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 3000)) %>% unique(), "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_plus_random_3000_100bp_tiles.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 5000)) %>% unique(), "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_plus_random_5000_100bp_tiles.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 50000)) %>% unique(), "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_plus_random_50000_100bp_tiles.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  
  
  
  #get annotation info
  mm10KG_txdb <- makeTxDbFromUCSC(genome="mm10", tablename="knownGene")
  bed_path <- file.path("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/", "mm10KnownGenes.bed")
  rtracklayer::export(asBED(mm10KG_txdb), bed_path)
  
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


no_prolong = function()
{
  meth_call_files_dir="C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/bismkark_meth_extractor_output_no_prolong"
  pipeline="bismarkCoverage"
  samp_ids=c("YoungYoung1", "YoungYoung2", "OldOld1", 
                        "Young1", "OldOld2", "Old1", "Young2",
                        "OldOld3", "Old2")
  treatments=c(1, 1, 1, 0, 1, 0, 0, 1, 0 )

  
  
  dir.create("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong")
  setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong")
  dir.create("figures")
  setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/figures")
  
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
  #
  # get hypo methylated bases
  dmrs_25p_hypo = getMethylDiff(tiles_raw_Cov10_unite_DMRs,difference=25,
                                qvalue=0.01,type="hypo")
  
  # visualize the distribution of hypo/hyper-methylated bases/regions per chromosome 
  diffMethPerChr(tiles_raw_Cov10_unite_DMRs,plot=TRUE,qvalue.cutoff=0.01,
                 meth.cutoff=25)
  
  #write methDiff files
  write.table(tiles_raw_Cov10_unite_DMRs, "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs.tsv",sep="\t")
  write.table(dmrs_25p_hyper, "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_25p_hyper.tsv",sep="\t")
  write.table(dmrs_25p_hypo, "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_25p_hypo.tsv",sep="\t")
  
  #write bed files (only chr start end)
  write.table(getData(dmrs_25p_hyper)[,1:3], "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_25p_hyper.bed",sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(getData(dmrs_25p_hypo)[,1:3], "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_25p_hypo.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(getData(tiles_raw_Cov10_unite)[,1:3], "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/all_100bp_tiles_united.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  #bg for graet
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 3000)) %>% unique(), "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_plus_random_3000_100bp_tiles.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 5000)) %>% unique(), "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_plus_random_5000_100bp_tiles.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  write.table(rbind(getData(dmrs_25p_hyper)[,1:3], getData(dmrs_25p_hypo)[,1:3], sample_n(getData(tiles_raw_Cov10_unite)[,1:3], 50000)) %>% unique(), "C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/noprolong/dmrs_plus_random_50000_100bp_tiles.bed",sep="\t",  row.names = FALSE , col.names = FALSE, quote = FALSE)
  
  
  
  #get annotation info
  mm10KG_txdb <- makeTxDbFromUCSC(genome="mm10", tablename="knownGene")
  bed_path <- file.path("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/", "mm10KnownGenes.bed")
  rtracklayer::export(asBED(mm10KG_txdb), bed_path)
  
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
