library(dmrseq)
library(dplyr)
library(BiocParallel)

register(MulticoreParam(6))


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

  for (dmrs in c("hypo_dmrs", "hyper_dmrs")) {
    bed_path <- paste0(output_path, "/", dmrs, ".bed")
    filtered_bed_path <-
      paste0(output_path, "/", dmrs, "_", min_meth_diff, "p.bed")
    ranked_path <- paste0(output_path, "/", dmrs, "_ranked.txt")

    # get dmrs as a data frame, sorted by standard bed order
    dmrs_df <- get(dmrs) %>% sort() %>% as.data.frame()

    # write.table(dmrs_df[, 1:3],
    #   file = bed_path, sep = "\t",
    #   row.names = FALSE, col.names = FALSE, quote = FALSE
    # )
    write_table(dmrs_df[, 1:3], bed_path)

    # write.table(dmrs_df[dmrs_df$meth_diff >= min_meth_diff, 1:3],
    #   file = filtered_bed_path, sep = "\t",
    #   row.names = FALSE, col.names = FALSE, quote = FALSE
    # )
    write_table(
            dmrs_df[dmrs_df$meth_diff >= min_meth_diff, 1:3], filtered_bed_path)

    # sort by descending width, descending diff, ascending qval
    dmrs_df <- dmrs_df %>%
      arrange(-width, -abs(meth_diff), qval)
  }
  write_table(dmrs_df, ranked_path)
  return(NULL)
}


# write all data in significant_regions as tsv and GRanges rds
write_full_dmr_data <- function(significant_regions, output_path) {
  # Function definition goes here
}



#' Find DMRs using dmrseq
#' @param samplesheet_path character, path to a samplesheet csv file
#'        see example_samplesheet.csv. the order of conditions determines the
#'        comparison direction, first vs second. I360 vs WT in the example.
#' @param output_path character, path to save bed files and other outputs to
#' @param min_qval numeric, minimum qval for a DMR to be considered significant
#' @param ... additional arguments passed to dmrseq
#' @return NULL
#' @details This function saves the identified DMRs as BED files.
#'          It also saves full DMRs info as a TSV file
#'          and an RDS file containg a GRanges object.
main <- function(samplesheet_path = "./samplesheet.csv", output_path = "./", min_qval = 0.05, ...) {
  # the covarite to test for differential methylation
  # harcoded because to change this, the function needs to be modified
  testCovariate <- "condition"


  # read samplesheet with files, conditions, replicates, (other covariates)
  samplesheet <- read.csv(samplesheet_path,
    header = TRUE, stringsAsFactors = FALSE
  )

  # build files vector ans sample description data frame from samplesheet
  samp_description_df <- data.frame(
    row.names = samplesheet$sample_name,
    condition = factor(
      samplesheet$condition,
      levels = samplesheet$condition %>% unique() %>% rev()
    ),
    replicate = samplesheet$replicate
  )

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

  # indeces of CpGs that have coverage in all samples
  loci_idx <- which(
    DelayedMatrixStats::rowSums2(getCoverage(bs, type = "Cov") == 0) == 0
  )

  # only work on chr1 for debugging
  loci_idx <- loci_idx[which(seqnames(bs)[loci_idx] == "chr1")] # debug

  # How many CpGs have coverage in all samples?
  sprintf(
    "covered CpGs: %d\ncovered in all samps: %d (%.2f%%)\n",
    nrow(bs), length(loci_idx), 100 * (length(loci_idx) / nrow(bs))
  ) %>% cat()

  # filter bs object to only keep CpGs with coverage in all samples
  bs <- bs[loci_idx, ]

  # run dmrseq
  regions <- dmrseq(bs, testCovariate, ...)

  significant_regions <- regions[regions$qval <= min_qval, ]

  n_hyper <- sum(significant_regions$stat > 0)
  n_hypo <- sum(significant_regions$stat < 0)
  sprintf(
    "number of hyper DMRs: %d\nnumber of hypo DMRs: %d\n",
    n_hyper, n_hypo
  ) %>% cat()


  # add % methylation difference column to significant_regions
  raw_diff <- meanDiff(bs, significant_regions, testCovariate)
  significant_regions$meth_diff <- 100 * raw_diff

  # write bed files of dmrs
  # hyper and hypo separately
  # standard bed sorted and sorted by width, diff, qval
  write_bed_files(significant_regions, output_path)

  # write all data in significant_regions as tsv and GRanges rds
  write_full_dmr_data(significant_regions, output_path)


  return(NULL)
}

setwd("/home/bengst/PycharmProjects/repo_for_reizel_lab/scripts_2024/dmrseq/")
main(samplesheet_path = "./example_samplesheet.csv", output_path = "./example_output", min_qval = 0.05)
