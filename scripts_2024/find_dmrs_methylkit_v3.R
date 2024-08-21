#!/usr/bin/env Rscript

#' find_dmrs_methylkit_v2_2023.R
#' Date: 30.07.2023
#' Author: Ben Steinberg
#' Purpose: find dmrs using methylkit (based on my script from 2022)

suppressMessages(library(argparser))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(methylKit))
# suppressMessages(library(GenomicFeatures)) # for getting annotation info
# suppressMessages(library(genomation)) #for annotating
# suppressMessages(library(rGREAT))
# library(R.utils) #for gunzip - not needed if not using the function get_KnownGene_from_ucsc()

#' Finding DMRs with methylKit
#'
#' @param meth_call_files_dir directory where the .cov files are (all will be used).
#' @param samp_ids character vector with the names of the samples (must match the order of the .cov files).
#' @param treatments  vector with the condition of each sample (0 or 1) the dmrs
#'                    are found as the difference between  1 - 0 groups (1 - treated , 0 - control).
#' @param pipeline name of the alignment pipeline, it can be either "amp", "bismark","bismarkCoverage", "bismarkCytosineReport" or a list (default:'amp'). See methylkit documentation for more details.
#' @param output_dir directory to save the results in.
#' @param meth_difference difference in percent for DMRs, default 25%.
#' @param genome refrence build e.g. mm9, mm10, hg38, etc.
#' @param base_cov minimum coverage per CpG, default 1.
#' @param tile_cov minimum coverage per tile, default 10.
#' @param tile_size tile size, default 100.
#' @param filt_hi_perc filter out bases with coverage above this percentile, default 99.9. Set to NULL to disable.
#' @param mc.cores number of cores to use for unite() and calculateDiffMeth(). must be set to 1 in Windows.
#'
#' @return
#' @export methdiff files , bed files, figures
#'
#' @examples
#' find_dmrs_main(meth_call_files_dir = "/meth_call_files_dir",
#'                samp_ids = c("LNCaP_1", "LNCaP_2", "LNCaP_3", "LNCaP_4", "LNCaP_5", "PrEC_1", "PrEC_2", "PrEC_3", "PrEC_4"),
#'                treatments = c(1, 1, 1, 1, 1, 0, 0, 0, 0),
#'                pipeline = "bismarkCoverage",
#'                output_dir = "/output_dir",
#'                meth_difference = 25,
#'                genome = "hg19",
#'                base_cov = 1,
#'                tile_cov = 10,
#'                tile_size = 100,
#'                filt_hi_perc = NULL,
#'                mc.cores = 10)
#'
find_dmrs_main <-
  function(meth_call_files_dir, samp_ids, treatments, pipeline, output_dir, meth_difference,
           genome, base_cov = 1, tile_cov = 10, tile_size = 100, filt_hi_perc = 99.9, mc.cores = 1)
  {
    set_up_directories(output_dir)
    print_args(genome, meth_call_files_dir, meth_difference,
               output_dir, pipeline, samp_ids, treatments, base_cov, tile_cov,
               tile_size, filt_hi_perc, mc.cores)

    tiles_raw_unite <- make_tiles(meth_call_files_dir, pipeline, samp_ids, treatments, genome,
                                  base_cov, tile_cov, tile_size, filt_hi_perc, mc.cores)
    plot_corelation_pca_hc(tiles_raw_unite)
    tiles_raw_unite_DMRs <- calculateDiffMeth(tiles_raw_unite, mc.cores = mc.cores)

    # get hyper and hypo methylated bases
    dmrs_hyper <- getMethylDiff(tiles_raw_unite_DMRs, difference = meth_difference, qvalue = 0.01, type = "hyper")
    dmrs_hypo <- getMethylDiff(tiles_raw_unite_DMRs, difference = meth_difference, qvalue = 0.01, type = "hypo")

    plot_meth_diff_per_chr(meth_difference, tiles_raw_unite_DMRs)

    # save objects. MethylKit objects cannot be loaded without methylKit package, objects are
    # therefore saved as GRanges objects by default, and optionally as MethylKit objects.
    info_str <- str_c(tile_size, "bpTiles_", meth_difference, "perc_", "basecov", base_cov, "_tilecov", tile_cov)
    save_methykit_objects(dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str, output_dir)
    save_granges_objects(dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str, output_dir)
    #TODO: add option to save as MethylKit objects

    #TODO: do we need meth scores in bed format? maybe bigWIg is better (can be done with rtracklayer)
    #write all tiles with meth scores (not only dmrs) that can be used for heatmaps and other applications
    write_meth_scores(tiles_raw_unite, str_c(output_dir, "/all_samps_", tile_size, "bp_tiles_meth_scores.bed"))

    write_bed_files(output_dir, dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str)
    write_bg_for_great(output_dir, dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str)


    #TODO: move annotation and gene ontology to separate script
    #NOTE: Must supply local KnownGenes.bed to get annotations.
    # Downloading KnownGenes.bed broken since 07-2023 (TxdDbFromUCSC is broken)
    # bed_path <- get_genes_info_old_broken_version(genome, known_genes_file, output_dir)
    # bed_path <- get_KnownGene_from_ucsc(genome, known_genes_file, output_dir) #TODO: figure out hot to convert KnownGenes.txt.gz to bed12
    # if (!is.null(known_genes_file)) {
    #   print("DEBUG: got known_genes_file")
    #   bed_path <- known_genes_file
    #   plot_dmr_gene_annotation(bed_path, dmrs_hyper, dmrs_hypo, meth_difference)
    # }
    # run_gene_ontology_analysis(dmrs_hypo, genome)
  }

parse_cli_args <- function() {
  p <- arg_parser("Find DMRs with methylKit")
  p <- add_argument(p, "--meth_call_files_dir", help = "directory where the .cov files are (all will be used)", short = "-m")
  p <- add_argument(p, "--samp_ids", help = "vector with the names of the samples separated by \"-\" (must match the order of the .cov files)", short = "-s")
  p <- add_argument(p, "--treatments", help = "vector with the condition of each sample (0 or 1) separated by \"-\" the dmrs are found as the difference between  1 - 0 groups (1 - treated , 0 - control)", short = "-t")
  p <- add_argument(p, "--pipeline", help = "name of the alignment pipeline, it can be either amp, bismark,bismarkCoverage, bismarkCytosineReport or a list. See methylkit documentation for more details.", short = "-p")
  p <- add_argument(p, "--output_dir", help = "directory to save the results in", short = "-o")
  p <- add_argument(p, "--genome", help = "mm9, mm10, hg38, etc.")
  p <- add_argument(p, "--known_genes_file", help = "annotation info in bed12 format e.g. mm10KnownGenes.bed.
                    For now must download manually from http://genome.ucsc.edu/cgi-bin/hgTables or use what I have
                    previously downloaded #TODO: if none is given will be downloaded - broken since 07-2023 (TxdDbFromUCSC is broken)")
  p <- add_argument(p, "--meth_difference", help = "difference in percent for DMRs, default 25%", default = 25)
  p <- add_argument(p, "--base_cov", help = "minimum coverage per CpG, default 1", default = 1)
  p <- add_argument(p, "--tile_cov", help = "minimum coverage per tile, default 10", default = 10)
  p <- add_argument(p, "--tile_size", help = "tile size, default 100", default = 100)
  p <- add_argument(p, "--filt_hi_perc", help = "filter out bases with coverage above this percentile, (#TODO: not sure this is needed for deduplicated WGBS)",
                    default = "99.9", type = "character") #type char to allow NULL
  p <- add_argument(p, "--mc.cores", help = "number of cores to use for unite() and calculateDiffMeth(). must be set to 1 in Windows", default = 1)
  argv <- parse_args(p)
  return(argv)
}

process_cli_args <- function(argv) {
  # to solve problem on condor multiple jobs will be using "-" instead of whitespace
  # to separate treatments and samp_ids
  # treatments <<- strsplit(argv$treatments, '-')[[1]] %>% as.numeric
  # samp_ids <<- strsplit(argv$samp_ids, '-')[[1]]

  # Not using condor anymore, so no need to split by "-". Just split by whitespace
  treatments <<- strsplit(argv$treatments, ' +')[[1]] %>% as.numeric
  samp_ids <<- strsplit(argv$samp_ids, ' +')[[1]]

  #allow list(...) as pipline input
  if (str_detect(argv$pipeline, "list"))
    argv$pipeline <<- eval(parse(text = argv$pipeline))

  if (str_detect(argv$filt_hi_perc, "NULL"))
    argv$filt_hi_perc <<- NULL
  else
    argv$filt_hi_perc <<- as.numeric(argv$filt_hi_perc)

  # not sure why normalizePath(argv$meth_call_files_dir) acts wiredly, if the argument is given straight to main it only
  # returns the relative path unless I print it first, loooks like a wiered bug.
  meth_call_files_dir <<- normalizePath(argv$meth_call_files_dir, winslash = "/", mustWork = TRUE)

  # dir must be created before running main because normalizePath()'s behavior is undefined if the dir dosn't exist
  create_dir(argv$output_dir)
  output_dir <<- normalizePath(argv$output_dir, winslash = "/", mustWork = TRUE)
}

read_meth_call_files <- function(meth_call_files_dir, pipeline_, samp_ids, treatments, genome, mincov = 10)
{
  meth_call_files <- list.files(path = meth_call_files_dir,
                                pattern = "*.cov.gz|*.cov",
                                full.names = TRUE)

  #try to use basename to get file names without full path, if that fails use str_split
  tryCatch(meth_call_files_no_fullpath <<- basename(meth_call_files),
           error = function(e) {
             print(e)
             split = str_split(meth_call_files, "/")
             meth_call_files_no_fullpath <<- lapply(split, function(x) x[length(x)]) %>% unlist()
           })

  sprintf("Make sure the samp_ids match the cov files order:\n") %>% cat()
  sprintf("cov file: %s\n samp id: %s\n\n", meth_call_files_no_fullpath, samp_ids) %>% cat()

  methyl_raw_list <- methRead(as.list(meth_call_files),
                              sample.id = as.list(samp_ids),
                              assembly = genome,
                              pipeline = pipeline_,
                              header = FALSE,
                              treatment = treatments,
                              context = "CpG",
                              mincov = mincov)
  return(methyl_raw_list)
}

filter_bases <- function(methyl_raw_list, lo.count = 10, hi.perc = 99.9)
{
  # It might be useful to filter samples based on coverage. Particularly, 
  # if our samples are suffering from PCR bias it would be useful to discard bases
  # with very high read coverage. Furthermore, we would also like to discard base
  # that have low read coverage, a high enough read coverage will increase the
  # power of the statistical tests. The code below filters a methylRawList and
  # discards bases that have coverage below "lo.count" and also discards the bases
  # that have more than 99.9th percentile of coverage in each sample.

  #NOTE: If we remove PCR duplicates e.g. with bismark deduplication, should hi.perc be set to 100 (or NULL)?

  filtered <- filterByCoverage(methyl_raw_list, lo.count = lo.count, lo.perc = NULL,
                               hi.count = NULL, hi.perc = hi.perc)
  return(filtered)
}

make_tiles <- function(meth_call_files_dir, pipeline, samp_ids, treatments, genome, base_cov = 1,
                       tile_cov = 10, tile_size = 100, filt_hi_perc = 99.9, mc.cores = 1)
{
  cat("Reading methylation files\n")
  methyl_raw_list <- read_meth_call_files(meth_call_files_dir, pipeline, samp_ids, treatments, genome, base_cov)

  #TODO: for deduplicated bismark change hi.prec to 100?
  cat("\nFiltering bases\n")
  methyl_raw_list <- filter_bases(methyl_raw_list, lo.count = base_cov, hi.perc = filt_hi_perc)
  # getMethylationStats(methyl_raw_list[[3]],plot=T, both.strands=FALSE)
  # getCoverageStats(methyl_raw_list[[3]],plot=T,both.strands=FALSE)

  cat("Making tiles\n")
  tiles_raw <- tileMethylCounts(methyl_raw_list, win.size = tile_size, step.size = tile_size,
                                mc.cores = mc.cores, cov.bases = tile_cov)

  # By default, unite function produces bases/regions covered in all samples.
  # That requirement can be relaxed using ???min.per.group??? option in unite function.
  cat("Uniting tiles from all samples into a single methylRaw object\n")
  tiles_raw_unite <- unite(tiles_raw, destrand = FALSE, mc.cores = mc.cores)

  return(tiles_raw_unite)
}

#' write all meth scores (not only dmrs) that can be used for heatmaps and other applications instead of
#' the tiles produced by the rrbs pipeline
write_meth_scores <- function(methylBase.obj, output_file)
{
  table <- percMethylation(methylBase.obj, rowids = TRUE)
  rownames(table) <- str_replace_all(rownames(table), '\\.', '\t')
  write.table(paste("#chr\tstart\tend", paste(colnames(table), collapse = "\t"), sep = "\t"),
              output_file, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  write.table(table / 100, output_file, append = TRUE, sep = "\t",
              row.names = TRUE, col.names = FALSE, quote = FALSE)
}

#' Plot correlation matrix, hierarchical clustering and PCA
plot_corelation_pca_hc <- function(tiles_raw_unite) {
  png(file = "correlation_matrix.png", width = 1000, height = 1000)
  correlation_matrix <- getCorrelation(tiles_raw_unite, plot = TRUE)
  dev.off()
  png(file = "hierarchical_clustering.png", width = 2000, height = 1000, res = 200)
  hc <- clusterSamples(tiles_raw_unite, dist = "correlation", method = "ward",
                       plot = TRUE)
  dev.off()
  png(file = "pca.png", width = 1500, height = 1000, res = 100)
  PCASamples(tiles_raw_unite)
  dev.off()
}

#TODO: fix or delete (used to work in 2022, broke in 2023)
get_genes_info_old_broken_version <- function(genome, known_genes_file, output_dir) {
  if (is.null(known_genes_file))
  {
    #TODO: on 26.07.2023 this is broken. knownGene has no track info.
    # Error in normArgTrack(track, trackids) (debug_makeTxDbFromUCSC.R#2): 'track' must be a single string or a list of strings

    #get annotation info
    #download  KnownGenes.bed file it wasn't given by user
    KG_txdb <- makeTxDbFromUCSC(genome = genome, tablename = "knownGene")
    bed_path <- file.path(output_dir, paste(genome, "KnownGenes.bed"))
    rtracklayer::export(asBED(KG_txdb), bed_path)
  }
  else
    bed_path <- known_genes_file
  return(bed_path)
}

#' Download KnownGenes.txt.gz file from UCSC and convert it to bed
#' Replaces get_genes_info_old_broken_version() that utalized makeTxDbFromUCSC
#' TODO: also broken - how to convert KnownGenes.txt.gz to bed12?
#' for now must download manually from http://genome.ucsc.edu/cgi-bin/hgTables
# get_KnownGene_from_ucsc <- function(genome, known_genes_file, output_dir) {
#   if (is.null(known_genes_file))
#   {
#     download.file(url = str_c("http://hgdownload.cse.ucsc.edu/goldenPath/", genome, "/database/knownGene.txt.gz"),
#                   destfile = str_c(output_dir, "/", genome, "KnownGenes.txt.gz"))
#     gunzip(str_c(output_dir, "/", genome, "KnownGenes.txt.gz"))
#     #TODO how the hell does one convert this to bed12?
#   }
#   else
#     bed_path <- known_genes_file
#   return(bed_path)
# }

#' Visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
#' @param meth_difference difference in percent for DMRs
#' @param tiles_raw_Cov10_unite_DMRs
plot_meth_diff_per_chr <- function(meth_difference, tiles_raw_unite_DMRs) {
  png("meth_diff_per_chr.png")
  diff_meth_per_chr <- diffMethPerChr(tiles_raw_unite_DMRs, plot = TRUE, qvalue.cutoff = 0.01,
                                      meth.cutoff = meth_difference)
  dev.off()
}

#' Save objects to disk
#' dmrs are saved as methylDiff objects
#' tiles methylation info is saved as a methylRaw object
save_methykit_objects <- function(dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str, output_dir) {
    saveRDS(dmrs_hyper, file = str_c(output_dir, "/dmrs_hyper_", info_str, ".rds.methylDiff"))
    saveRDS(dmrs_hypo, file = str_c(output_dir, "/dmrs_hypo_", info_str, ".rds.methylDiff"))
    saveRDS(tiles_raw_unite, file = str_c(output_dir, "/raw_unite_", info_str, ".rds.methylRaw"))
  }

save_granges_objects <- function(dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str, output_dir) {
  saveRDS(as(dmrs_hyper, "GRanges"), file = str_c(output_dir, "/dmrs_hyper_", info_str, ".rds.GRanges"))
  saveRDS(as(dmrs_hypo, "GRanges"), file = str_c(output_dir, "/dmrs_hypo_", info_str, ".rds.GRanges"))
  saveRDS(as(tiles_raw_unite, "GRanges"), file = str_c(output_dir, "/raw_unite_", info_str, ".rds.GRanges"))
}

#' Write bed files (only chr start end)
#' tiles_raw_unite - all tiles covered by all samples
write_bed_files <- function(output_dir, dmrs_hyper, dmrs_hypo, tiles_raw_unite, info_str) {
  write.table(getData(dmrs_hyper)[, 1:3], str_c(output_dir, "/dmrs_hyper_", info_str, ".bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(getData(dmrs_hypo)[, 1:3], str_c(output_dir, "/dmrs_hypo_", info_str, ".bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(getData(tiles_raw_unite)[, 1:3], str_c(output_dir, "/raw_unite_", info_str, ".bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' Write backgronud files for GREAT (gene ontology analysis) web GUI
write_bg_for_great <- function(output_dir, dmrs_hyper, dmrs_hypo, tiles_raw_unite,
                               info_str, n_samples = c(3000, 5000, 50000)) {
  if (nrow(dmrs_hyper) == 0) {
    cat("No hypermethylated DMRs found - background BED files for GREAT not written\n")
    return(1)
  }
  if (nrow(dmrs_hypo) == 0) {
    cat("No hypomethylated DMRs found - background BED files for GREAT not written\n")
    return(1)
  }
  if (nrow(tiles_raw_unite) == 0) {
    cat("No tiles with coverage in all samples - background BED files for GREAT not written\n")
    return(1)
  }

  # write.table(rbind(getData(dmrs_hyper)[, 1:3], getData(dmrs_hypo)[, 1:3],
  #                   sample_n(getData(tiles_raw_unite)[, 1:3], 3000)) %>% unique(),
  #             str_c(output_dir, "/dmrs_plus_random_3000_", info_str, ".bed"),
  #             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.table(rbind(getData(dmrs_hyper)[, 1:3], getData(dmrs_hypo)[, 1:3],
  #                   sample_n(getData(tiles_raw_unite)[, 1:3], 5000)) %>% unique(),
  #             str_c(output_dir, "/dmrs_plus_random_5000_", info_str, ".bed"),
  #             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.table(rbind(getData(dmrs_hyper)[, 1:3], getData(dmrs_hypo)[, 1:3],
  #                   sample_n(getData(tiles_raw_unite)[, 1:3], 50000)) %>% unique(),
  #             str_c(output_dir, "/dmrs_plus_random_50000_", info_str, ".bed"),
  #             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # replace above block with loop
  for (n in n_samples) {
    if (nrow(tiles_raw_unite) < n) {
      cat(str_c("Not enough regions to randomly sample n=", n
        , " from  - background BED file for GREAT not written\n"))
      next
    }
    write.table(rbind(getData(dmrs_hyper)[, 1:3], getData(dmrs_hypo)[, 1:3],
                      sample_n(getData(tiles_raw_unite)[, 1:3], n)) %>% unique(),
                str_c(output_dir, "/dmrs_plus_random_", n, "_", info_str, ".bed"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

#' Plot pie chart with relative distribution of exons, introns, promoters, intergenic regions ..
plot_dmr_gene_annotation <- function(bed_path, dmrs_hyper, dmrs_hypo, meth_difference) {
  gene.obj <- readTranscriptFeatures(bed_path)
  dmrs_hyper_annotation <- annotateWithGeneParts(as(dmrs_hyper, "GRanges"), gene.obj)
  dmrs_hypo_annotation <- annotateWithGeneParts(as(dmrs_hypo, "GRanges"), gene.obj)
  png(file = "hypo_annotation.png", width = 1000, height = 1000, res = 150)
  plotTargetAnnotation(dmrs_hypo_annotation, precedence = TRUE, main = str_c("DMRs ", meth_difference, "% hypo annotation"))
  dev.off()
  png(file = "hyper_annotation.png", width = 1000, height = 1000, res = 150)
  plotTargetAnnotation(dmrs_hyper_annotation, precedence = TRUE, main = str_c("DMRs ", meth_difference, "% hyper annotation"))
  dev.off()
}

#'Gene Ontology analysis via GREAT
run_gene_ontology_analysis <- function(dmrs_hypo, genome) {
  #TODO: ask tzachi about great params, use bg?
  great_job <- submitGreatJob(as(dmrs_hypo, "GRanges"), species = genome)
  tb <- getEnrichmentTables(great_job, download_by = "tsv")
  png(file = "RegionGeneAssociationGraphs.png", width = 10000, height = 2500, res = 1000)
  res = plotRegionGeneAssociationGraphs(great_job)
  dev.off()
}

print_args <- function(genome, meth_call_files_dir, meth_difference, output_dir, pipeline, samp_ids,
                       treatments, base_cov, tile_cov, tile_size, filt_hi_perc, mc.cores) {
  sprintf("\n\nGiven Arguments:\n") %>% cat()
  sprintf("pwd: %s\n", getwd()) %>% cat()
  sprintf("meth_call_files_dir: %s\n", meth_call_files_dir) %>% cat()
  sprintf("output_dir: %s\n", output_dir) %>% cat()
  sprintf("samp_ids: %s\n", str_c(samp_ids, collapse = ", ")) %>% cat()
  sprintf("treatments: %s\n", str_c(treatments, collapse = ", ")) %>% cat()
  sprintf("pipeline: %s\n", pipeline) %>% cat()
  sprintf("meth_difference: %s\n", meth_difference) %>% cat()
  sprintf("base_cov: %s\n", base_cov) %>% cat()
  sprintf("tile_cov: %s\n", tile_cov) %>% cat()
  sprintf("tile_size: %s\n", tile_size) %>% cat()
  sprintf("genome: %s\n", genome) %>% cat()
  sprintf("filt_hi_perc: %s\n", filt_hi_perc) %>% cat()
  sprintf("mc.cores: %s\n\n\n", mc.cores) %>% cat()
}

create_dir <- function(output_dir) {
  if (!dir.exists(output_dir))
    dir.create(output_dir)
}


#' make sure output_dir exists, create figures dir inside it, and setwd to figures
set_up_directories <- function(output_dir) {
  if (!exists(output_dir)) #if using main() directly
    create_dir(output_dir)
  setwd(output_dir)
  if (!dir.exists("figures"))
    dir.create("figures")
  setwd("./figures")
}


#TODO: find differentialy methylated cytosines or CpGs, possibly group them by reigon to find DMRs

##########
## main ##
##########

# Create a parser
argv <- parse_cli_args()

#print cli args
cat("\n\nRunning:\n", commandArgs(trailingOnly = FALSE), "\n")
process_cli_args(argv)
find_dmrs_main(meth_call_files_dir, samp_ids, treatments, argv$pipeline,
               output_dir, as.numeric(argv$meth_difference), argv$genome,
               as.numeric(argv$base_cov), as.numeric(argv$tile_cov), argv$tile_size, argv$filt_hi_perc,
               as.numeric(argv$mc.cores)
)

#use main manually
# find_dmrs_main(meth_call_files_dir = "/meth_call_files_dir",
#                samp_ids = c("LNCaP_1", "LNCaP_2", "LNCaP_3", "LNCaP_4", "LNCaP_5", "PrEC_1", "PrEC_2", "PrEC_3", "PrEC_4"),
#                treatments = c(1, 1, 1, 1, 1, 0, 0, 0, 0),
#                pipeline = "bismarkCoverage",
#                output_dir = "/output_dir",
#                known_genes_file = "",
#                meth_difference = 25,
#                genome = "hg19",
#                base_cov = 10,
#                tile_cov = 10,
#                tile_size = 100,
#                filt_hi_perc = NULL,
#                mc.cores = 10
# )

# meth_call_files_dir = "/meth_call_files_dir"
# samp_ids = c("LNCaP_1", "LNCaP_2", "LNCaP_3", "LNCaP_4", "LNCaP_5", "PrEC_1", "PrEC_2", "PrEC_3", "PrEC_4")
# treatments = c(1, 1, 1, 1, 1, 0, 0, 0, 0)
# pipeline = "bismarkCoverage"
# output_dir = "/output_dir"
# known_genes_file = ""
# meth_difference = 25
# genome = "hg19"
# base_cov = 10
# tile_cov = 10
# tile_size = 100
# filt_hi_perc = 0
# mc.cores = 1