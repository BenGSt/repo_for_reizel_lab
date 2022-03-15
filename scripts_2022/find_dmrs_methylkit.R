library(methylKit)


install_methylkit = function()
{
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("methylKit")
  
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


meth_call_files_dir="C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/bismkark_meth_extractor_output"
pipeline="bismarkCoverage"
samp_ids=c("YoungYoung1", "YoungYoung2", "OldOld1", 
           "Young1", "OldOld2", "Old1", "Young2",
           "OldOld3", "Old2")
treatments=c(1,1,1,0,1,0,0,1,0)

methyl_raw_list = read_meth_call_files(
                            meth_call_files_dir,
                            pipeline,
                            samp_ids,
                            treatments)

methyl_raw_list = filter_bases(methyl_raw_list)
getMethylationStats(methyl_raw_list[[3]],plot=T, both.strands=FALSE)
getCoverageStats(methyl_raw_list[[3]],plot=T,both.strands=FALSE)

tiles_raw = tileMethylCounts(methyl_raw_list, win.size=100, step.size=100)
tiles_raw_Cov10_unite=unite(tiles_raw, destrand=FALSE)

tiles_raw_Cov10_unite_DMRs=calculateDiffMeth(tiles_raw_Cov10_unite)
write.table(tiles_raw_Cov10_unite_DMRs,"C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/dmrs_no_YYProlong.tsv",sep="\t")

