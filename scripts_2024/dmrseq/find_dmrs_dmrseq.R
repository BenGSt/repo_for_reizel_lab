write_table <- function(to_write, path) {
  write.table(to_write,
    file = path, sep = "\t",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

#' Write BED files of DMRs
#'
#' Writes BED format files for differentially methylated regions (DMRs),
#' separating hypermethylated and hypomethylated regions. Creates multiple
#' sorted versions of the files.
#'
#' @param significant_regions GRanges, dmrseq output + meth_diff % column
#' @param output_path character, directory to save the BED files to
#' @param min_meth_diff numeric, minimum absolute methylation difference
#'        [%] required for a region to be included in the filtered BED files
#' @details Creates separate BED files for hypermethylated (meth_diff > 0) and
#'          hypomethylated (meth_diff < 0) regions. For each type, generates:
#'          1. Standard sorted BED
#'             (alphabetically by chromosome and numerically by position)
#'          2. BED sorted by region width, methylation difference, and q-value
#' @return No return value, writes BED files to disk
write_bed_files <- function(
    significant_regions, output_path, min_meth_diff = 25) {
  if (!dir.exists(output_path)) {
    message("Creating output directory: ", output_path)
    dir.create(output_path, recursive = TRUE)
  }

  hypo_dmrs <- significant_regions[significant_regions$meth_diff < 0, ]
  hyper_dmrs <- significant_regions[significant_regions$meth_diff > 0, ]

  # print # of significant dmrs
  print_num_dmrs(nrow(hyper_dmrs), nrow(hypo_dmrs), significant = TRUE)

  for (dmrs in c("hypo_dmrs", "hyper_dmrs")) {
    bed_path <- paste0(output_path, "/", dmrs, ".bed")
    filtered_bed_path <-
      paste0(output_path, "/", dmrs, "_", min_meth_diff, "p.bed")
    ranked_path <- paste0(output_path, "/", dmrs, "_ranked.txt")

    # get dmrs as a data frame, sorted by standard bed order
    dmrs_df <- get(dmrs) %>%
      sort() %>%
      as.data.frame()

    write_table(dmrs_df[, 1:3], bed_path)

    write_table(
      dmrs_df[dmrs_df$meth_diff >= min_meth_diff, 1:3], filtered_bed_path
    )

    # sort by descending width, descending diff, ascending qval
    dmrs_df <- dmrs_df %>%
      arrange(-width, -abs(meth_diff), qval)
    write_table(dmrs_df, ranked_path)
  }

  return(NULL)
}


#' Write Full DMR Data to Files
#'
#' This function writes the significant regions data to a tsv file
#'
#' @param significant_regions A GRanges object containing the
#'        significant regions data.
#' @param output_path A string specifying the directory where the output files
#'        will be saved.
#'
#' @return NULL, writes files to disk.
write_full_dmr_data <- function(significant_regions, output_path, save_rds = FALSE) {
  full_data_tsv <- file.path(output_path, "full_dmr_data.tsv")
  write_table(significant_regions, full_data_tsv)
}


print_num_dmrs <- function(n_hyper, n_hypo, significant = FALSE) {
  prefix <- if (significant) "significant " else ""
  sprintf(
    "number of %shyper DMRs: %d\nnumber of %shypo DMRs: %d\n",
    prefix, n_hyper, prefix, n_hypo
  ) %>% cat()
  cat("\n")
}


#' Inform on how many CpGs have coverage in all samples
print_cpg_cov_info <- function(bs, loci_idx) {
  sprintf(
    "covered CpGs: %d\ncovered in all samps: %d (%.2f%%)\n",
    nrow(bs), length(loci_idx), 100 * (length(loci_idx) / nrow(bs))
  ) %>% cat()
  cat("\n")
}


#' Make a sample description data frame from a samplesheet
make_samp_description_df <- function(samplesheet) {
  samp_description_df <- data.frame(
    row.names = samplesheet$sample_name,
    condition = factor(
      samplesheet$condition,
      levels = samplesheet$condition %>% unique() %>% rev()
    ),
    replicate = samplesheet$replicate
  )
  return(samp_description_df)
}


#' Find DMRs using dmrseq
#' @param samplesheet_path character, path to a samplesheet csv file
#'        see example_samplesheet.csv. the order of conditions determines the
#'        comparison direction, first vs second. I360 vs WT in the example.
#' @param output_path character, path to save bed files and other outputs to
#' @param min_qval numeric, minimum qval for a DMR to be considered significant
#' @param min_meth_diff numeric, minimum absolute methylation difference [%] to write to 'filtered' bed file
#' @param save_full_dmr_data logical, whether to save the full DMR data icluding statistics as a TSV file
#' @param save_rds logical, whether to save the resulting objects as RDS files (can be loaded individually in R)
#' @param save_rdata logical, whether to save the R environment (can be loaded with load())
#' @param ... additional arguments passed to dmrseq
#' @return NULL
#' @details This function saves the identified DMRs as BED files.
#'          It also saves full DMRs info as a TSV file
#'          and an RDS file containg a GRanges object.
main <- function(
    samplesheet_path = "./samplesheet.csv",
    output_path = "./",
    min_qval = 0.05,
    min_meth_diff = 25,
    save_full_dmr_data = TRUE,
    save_rds = FALSE,
    save_rdata = FALSE,
    ...) {
  # the covarite to test for differential methylation:
  # harcoded because to change this the function needs to be modified
  testCovariate <- "condition"

  # read samplesheet with files, conditions, replicates, (other covariates)
  samplesheet <- read.csv(samplesheet_path,
    header = TRUE, stringsAsFactors = FALSE
  )

  # build files vector ans sample description data frame from samplesheet
  samp_description_df <- make_samp_description_df(samplesheet)
  files <- samplesheet$file_path

  # create a BSseq object, assume bismark cov files
  bs <- read.bismark(
    files = files,
    colData = samp_description_df,
    rmZeroCov = TRUE,
    strandCollapse = TRUE,
    verbose = TRUE
  )
  # Note: I don't think strandCollapse works with bismark cov
  # (can combine with bedgraph or bam + methylkit to get strand info and test)

  cat("bs object created\npData:\n")
  print(pData(bs))
  cat("\nbs:\n")
  print(bs)
  cat("\n")

  # indeces of CpGs that have coverage in all samples
  loci_idx <- which(
    DelayedMatrixStats::rowSums2(getCoverage(bs, type = "Cov") == 0) == 0
  )

  # only work on one chromosome for debugging
  loci_idx <- loci_idx[which(seqnames(bs)[loci_idx] == "chr3")] # debug

  # How many CpGs have coverage in all samples?
  print_cpg_cov_info(bs, loci_idx)

  # filter bs object to only keep CpGs with coverage in all samples
  bs <- bs[loci_idx, ]

  # run dmrseq
  regions <- dmrseq(bs, testCovariate, ...)
  significant_regions <- regions[regions$qval <= min_qval, ]

  n_hyper <- sum(significant_regions$stat > 0)
  n_hypo <- sum(significant_regions$stat < 0)
  print_num_dmrs(n_hyper, n_hypo)

  if (n_hyper + n_hypo > 0) {
    # add % methylation difference column to significant_regions
    raw_diff <- meanDiff(bs, significant_regions, testCovariate)
    significant_regions$meth_diff <- 100 * raw_diff

    # write bed files of dmrs, hyper and hypo separately
    # standard bed sorted and sorted by width, diff, qval
    write_bed_files(significant_regions, output_path, min_meth_diff)

    # write all data in significant_regions as tsv
    if (save_full_dmr_data) {
      write_full_dmr_data(significant_regions, output_path)
    }
  }

  if (save_rds) {
    if (n_hyper + n_hypo > 0) {
      full_data_rds <- file.path(output_path, "full_dmr_data.rds")
      saveRDS(significant_regions, full_data_rds)
    }
    saveRDS(regions, file.path(output_path, "dmrseq_regions.rds"))
    saveRDS(bs, file.path(output_path, "bs.rds"))
  }

  if (save_rdata) {
    save.image(file.path(output_path, "dmrseq_results.RData"))
  }

  return()
}



help <- function() {
  cat("Usage: Rscript find_dmrs_dmrseq.R -samplesheet <samplesheet_path> -output_dir <output_path> -min_qval <min_qval> [flags] [dmrseq_args]\n")
  cat("Arguments:\n")
  cat("-help: print this help message\n")
  cat("-ncpus: number of cores to use (default: 6)\n")
  cat("-samplesheet: path to a samplesheet csv file. (default:./samplesheet.csv)\n")
  cat("-output_dir: path to save bed files and other outputs to (defualt: ./dmrseq_output)\n")
  cat("-min_qval: minimum qval for a DMR to be considered significant (defualt: 0.05)\n")
  cat("-min_meth_diff: minimum absolute methylation difference [%] to write to 'filtered' bed file (default: 25)\n")
  cat("-save_full_dmr_data: flag, save the full DMR data including statistics as a TSV file\n")
  cat("-save_rds: flag, save the resulting objects as RDS files (can be loaded individually in R)\n")
  cat("-save_rdata: flag, save the R environment (can be loaded with load())\n")
  cat("[dmrseq_args]: additional arguments passed to dmrseq\n")
}


# set default params
args_list <- list()
ncpus <- 6
args_list[["samplesheet_path"]] <- "./samplesheet.csv"
args_list[["output_path"]] <- "./dmrseq_output"
args_list[["min_qval"]] <- 0.05
args_list[["min_meth_diff"]] <- 25
args_list[["save_full_dmr_data"]] <- FALSE
args_list[["save_rds"]] <- FALSE
args_list[["save_rdata"]] <- FALSE



# get params from cli
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  help()
  quit(status = 1)
}

while (length(args) > 0) {
  if (args[1] == "-help") {
    help()
    quit(status = 0)
  } else if (args[1] == "-ncpus") {
    ncpus <- as.numeric(args[2])
    args <- args[-(1:2)]
  } else if (args[1] == "-samplesheet") {
    args_list[["samplesheet_path"]] <- args[2]
    args <- args[-(1:2)]
  } else if (args[1] == "-output_dir") {
    args_list[["output_path"]] <- args[2]
    args <- args[-(1:2)]
  } else if (args[1] == "-min_qval") {
    args_list[["min_qval"]] <- as.numeric(args[2])
    args <- args[-(1:2)]
  } else if (args[1] == "-min_meth_diff") {
    args_list[["min_meth_diff"]] <- as.numeric(args[2])
    args <- args[-(1:2)]
  } else if (args[1] == "-save_full_dmr_data") {
    args_list[["save_full_dmr_data"]] <- TRUE
    args <- args[-1]
  } else if (args[1] == "-save_rds") {
    args_list[["save_rds"]] <- TRUE
    args <- args[-1]
  } else if (args[1] == "-save_rdata") {
    args_list[["save_rdata"]] <- TRUE
    args <- args[-1]
  } else if (args[1] %in% names(formals(dmrseq::dmrseq))) {
    type <- typeof(formals(dmrseq::dmrseq)[[args[1]]])
    args_list[[args[1]]] <- as(args[2], type)
    args <- args[-(1:2)]
  } else {
    stop("Unknown argument: ", args[1])
    help()
    quit(status = 1)
  }
}

library(dmrseq)
library(dplyr)
library(BiocParallel)

# setwd("/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/dmrseq/")
register(MulticoreParam(ncpus))

cat("Running with the following parameters:\n")
# print each arg in args_list nicely
cat("ncpus: ", ncpus, "\n")
for (arg in names(args_list)) {
  cat(arg, ": ", args_list[[arg]], "\n")
}


returned <- do.call(main, args_list)
# main(samplesheet_path = "./example_samplesheet.csv", output_path = "./example_output", min_qval = 0.05)
