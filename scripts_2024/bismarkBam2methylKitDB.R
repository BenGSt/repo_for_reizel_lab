#!/usr/bin/Rscript

# File: bismarkBam2methylKitDB.R
# Description: Script to process Bismark BAM files and convert to methylKit database
# Author: Ben Steinberg
# Date: 2024-09-18
# Version: 1.0

library(argparser)
# Define command line arguments
parser <- arg_parser("Process Bismark BAM files and convert to methylKit database")
parser <- add_argument(parser, "--bam_file", help = "Bismark BAM file to process", type = "character")
parser <- add_argument(parser, "--sample_id", help = "Sample ID", type = "character")
parser <- add_argument(parser, "--assembly", help = "Genome assembly (e.g., hg19)", type = "character")
parser <- add_argument(parser, "--save_folder", help = "Folder to save output files", type = "character")
parser <- add_argument(parser, "--write_cov_file", help = "Should a coverage file be saved",type = "logical", default = FALSE)
parser <- add_argument(parser, "--save_context", help = "Context to save to coverage file. One of: 'CpG','CHG','CHH'", type = "character", default = "CpG")
parser <- add_argument(parser, "--read_context", help = "Context to save to methylRawDB bgz file. One of: 'CpG','CHG','CHH' ,'none'", type = "character", default = "CpG")
parser <- add_argument(parser, "--nolap", help = "if set to TRUE and the BAM file has paired-end reads, the one read of the overlapping paired-end read pair will be ignored for methylation calling.", type = "logical", default = FALSE)
parser <- add_argument(parser, "--mincov", help = "Minimum coverage required to include a C in the output", type = "numeric", default = 1)
parser <- add_argument(parser, "--minqual", help = "Minimum quality required to include a C in the output", type = "numeric", default = 20)

argv <- parse_args(parser)

#check if arguments without default values are provided
if (is.na(argv$bam_file) || is.na(argv$sample_id) || is.na(argv$assembly) || is.na(argv$save_folder)) {
  stop("Please provide all required arguments.
  -b, --bam_file      Bismark BAM file to process
  -s, --sample_id     Sample ID
  -a, --assembly      Genome assembly (e.g., hg19)
  --save_folder       Folder to save output files
")
}

#print arguments nicely
cat("Arguments provided:
  Bismark BAM file: ", argv$bam_file, "
  Sample ID: ", argv$sample_id, "
  Genome assembly: ", argv$assembly, "
  Save folder: ", argv$save_folder, "
  Save coverage file: ", argv$write_cov_file, "
  Save context: ", argv$save_context, "
  Read context: ", argv$read_context, "
  Nolap: ", argv$nolap, "
  Minimum coverage: ", argv$mincov, "
  Minimum quality: ", argv$minqual, "\n")

# If we don't want to save a coverage file, we will cd to save_folder
# and methylKit will save the DB files there, in a directory
# named after this scheme: "methylDB <Date> <3randomlettersornumbers>".
# The files can be moved later.
# Annoying! but easier than rewriting methylKit::ProcessBismarkAln()
if (!argv$write_cov_file) {
  if (!dir.exists(argv$save_folder)) {
    dir.create(argv$save_folder)
  }
  #get full path to bam_file before changing directory
  argv$bam_file <- normalizePath(argv$bam_file)

  setwd(argv$save_folder)
  argv$save_folder <- NULL
}

library(methylKit)
# Process single Bismark BAM file and convert to coverage file and methylKit database
# Note: A DB can't be saved without saving the coverage file (save.context can't be "none"),
#       and reading the data to memory (the data defined in read.context is saved to DB).

# Process single Bismark BAM file and convert to coverage file and methylKit database
processBismarkAln(
  location = argv$bam_file,
  sample.id = argv$sample_id,
  assembly = argv$assembly,
  save.folder = argv$save_folder,
  save.context = argv$save_context,
  read.context = argv$read_context,
  nolap = argv$nolap,
  mincov = argv$mincov,
  minqual = argv$minqual,
  phred64 = FALSE,
  treatment = NULL,
  save.db = TRUE,
  verbose = 1
)

# #test
# bam_file <- "/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/test_bismarkBam2MethylKit/HFF_TLarge_ctrl_virusRatio.1MReads.sorted.bam"
# file.exists(bam_file)

# # Process single Bismark BAM file and convert to coverage file and methylKit database
# # Note: A DB can't be saved without saving the coverage file (save.context can't be "none"),
# #       and reading the data to memory (the data defined in read.context is saved to DB).
# setwd("/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/test_bismarkBam2MethylKit/")
# save_folder <- "/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/test_bismarkBam2MethylKit/DB/"
# if (!dir.exists(save_folder)) {
#   dir.create(save_folder)
# }
# setwd(save_folder)

# processBismarkAln(
#   location = bam_file,
#   sample.id = "HFF_TLarge_ctrl_virusRatio_1",
#   assembly = "hg19",
#   # save.folder = "/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/test_bismarkBam2MethylKit/DB/",
#   save.context = "CpG",
#   read.context = "CpG",
#   nolap = FALSE,
#   mincov = 1,
#   minqual = 20,
#   phred64 = FALSE,
#   treatment = NULL,
#   save.db = TRUE,
#   verbose = 1
# )



# How to use the produced files
example <- function(){
  # read the coverage file
  methyl_raw <- methRead(
    location = "/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/test_bismarkBam2MethylKit/DB/HFF_TLarge_ctrl_virusRatio_1_CpG.txt",
    sample.id = "HFF_TLarge_ctrl",
    assembly = "hg19",
    header = TRUE,
    pipeline = "bismark"
  )

  # read the database
  bgz_file <- "/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/test_bismarkBam2MethylKit/DB/HFF_TLarge_ctrl_virusRatio_1_cpg.txt.bgz"
  files <- c(bgz_file, bgz_file)
  list_meth_raws <- lapply(files, function(x) readMethylDB(x)[])

  meth_raw_list <- methylRawList(
    list_meth_raws,
    treatment = c(1, 0)
  )

  # try CpG merging from both starnds
  united <- unite(meth_raw_list, destrand = FALSE)
  united_destranded <- unite(meth_raw_list, destrand = TRUE)
  library(dplyr)

  # print the number of CpGs with coverage1 > 10
  getData(united) %>%
    filter(coverage1 > 10) %>%
    nrow()

  getData(united_destranded) %>%
    filter(coverage1 > 10) %>%
    nrow()


  brks <- seq(50, 100, 10)
  getData(united) %>%
    filter(coverage1 > 50, coverage1 < 100) %>%
    select(coverage1) %>%
    unlist() %>%
    hist(., breaks = brks, xaxs = "i", main = "United", xlab = "Coverage", ylab = "Frequency",
    col = rgb(0, 0, 1, 0.5))

  getData(united_destranded) %>%
    filter(coverage1 > 50, coverage1 < 100) %>%
    select(coverage1) %>%
    unlist() %>%
    hist(., breaks = brks, main = "United Destranded", xlab = "Coverage", ylab = "Frequency",
    add = TRUE, col = rgb(1, 0, 0, 0.5))
}