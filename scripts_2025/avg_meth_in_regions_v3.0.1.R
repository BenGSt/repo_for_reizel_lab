#!/usr/bin/env Rscript

#' --- version 3.0 chnages:
#' - Use BiocParallel for parallelization instead of parallel package
#' - increase memory efficiency by avoiding repeated object conversion in regionCounts
#' So far this works well. comitting here before making more changes.
#' #TODO: -- avoid accumulating counts in memory and then calculating avg methylation per tile,
#'           instead calculate avg methylation per tile on the fly and only keep the result

#' @title Average methylation by position in genomic regions
#' @description calculate the average methylation by position over a set genomic regions
#' by finding the methylation level in sliding window tiles for each region,
#' then taking the average of each tile over all regions.
#' @author Ben Steinberg, Pleleg Shalev
#' @date 6.9.2024
# TODO: tested on ido goldstein's data - set header = TRUE in methRead for methyldackel files - see TODOs in body

library(rtracklayer)
library(data.table)
library(dplyr)
library(BiocParallel)
library(argparser)
library(IRanges)

library(GenomeInfoDb)
library(GenomicRanges)
library(methylKit)

print_mem_usage <- function(stage) {
  cat(sprintf("DEBUG MEMORY [%s]:\n", stage))
  print(gc())
}

################################
# new/altered methylKit methods#
################################

# methylKit does not define a method "percMethylation" for "methylRaw"
# objects, only for "methylBase" objects. These need to contain more than one sample
# (experiment). Adding the method for the use case of region counts on one sample (used here for
# parallelization (so this script can be run on many samples in parallel)
setMethod(
  "percMethylation", "methylRaw",
  function(methylBase.obj, rowids = FALSE, save.txt = FALSE, chunk.size = 1e6) {
    meth <- 100 * methylBase.obj$numCs / (methylBase.obj$numCs + methylBase.obj$numTs)
    meth.mat <- as.matrix(meth)
    # names(meth.mat) <- methylBase.obj@sample.id
    # rownames(meth.mat) <- NULL
    # if (rowids) {
    #   rownames(meth.mat) <- as.character(paste(x[, 1], x[, 2], x[, 3], sep = "."))
    # }
    return(as.matrix(meth.mat))
  }
)

# Ben 1.9.2024: I need regionCounts to return all regions, even if they have no coverage
# I will modify the method:

# GETs regional counts for given GRanges object
# RETURNS a new methylRaw object
# @param object a \code{methylRaw} object
# @param regions a GRanges object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRaw,GRanges-method
setMethod(
  "regionCounts", signature(object = "methylRaw", regions = "GRanges"),
  function(object, regions, cov.bases, strand.aware, save.db = FALSE, ...) {
    # require(GenomicRanges)
    # sort regions
    regions <- sortSeqlevels(regions)
    regions <- sort(regions, ignore.strand = TRUE)
    
    # Check for pre-converted GRanges in ...
    args <- list(...)
    if ("g.meth" %in% names(args)) {
      g.meth <- args$g.meth
    } else {
      # overlap object with regions
      # convert object to GRanges
      g.meth <- as(object, "GRanges")
      if (!strand.aware) {
        strand(g.meth) <- "*"
      }
    }

    mat <- IRanges::as.matrix(findOverlaps(regions, g.meth))
    
    # find the regions that have no coverage
    regions_no_hits_indeces <- setdiff(seq_along(regions), unique(mat[, 1]))

    # create a temporary data.table row ids from regions and counts from object
    coverage <- numCs <- numTs <- id <- covered <- NULL
    df <- data.frame(id = mat[, 1], getData(object)[mat[, 2], c(5, 6, 7)])
    dt <- data.table(df)

    # use data.table to sum up counts per region
    sum.dt <- dt[, list(
      coverage = sum(coverage),
      numCs = sum(numCs),
      numTs = sum(numTs), covered = length(numTs)
    ), by = id]
    sum.dt <- sum.dt[sum.dt$covered >= cov.bases, ]

    # add regions with no coverage to sum.dt
    if (length(regions_no_hits_indeces) > 0) {
      sum.dt <- rbind(sum.dt, data.table(
        id = regions_no_hits_indeces,
        coverage = 0, numCs = NA, numTs = NA, covered = 0
      ))
      # sort by id to keep the order of the regions
      setkey(sum.dt, id)
    }

    temp.df <- as.data.frame(regions) # get regions to a dataframe

    # create id string for the new object to be returned
    # ids have to be unique and we can not assume GRanges objects will
    # have a name attribute
    if ("name" %in% names(temp.df)) {
      new.ids <- paste(temp.df[sum.dt$id, "seqnames"], temp.df[sum.dt$id, "start"],
        temp.df[sum.dt$id, "end"], temp.df[sum.dt$id, "name"],
        sep = "."
      )
    } else {
      new.ids <- paste(temp.df[sum.dt$id, "seqnames"], temp.df[sum.dt$id, "start"],
        temp.df[sum.dt$id, "end"],
        sep = "."
      )
    }

    # create a new methylRaw object to return
    new.data <- data.frame( # id      =new.ids,
      chr = temp.df[sum.dt$id, "seqnames"],
      start = temp.df[sum.dt$id, "start"],
      end = temp.df[sum.dt$id, "end"],
      strand = temp.df[sum.dt$id, "strand"],
      coverage = sum.dt$coverage,
      numCs = sum.dt$numCs,
      numTs = sum.dt$numTs
    )

    if (!save.db) {
      new("methylRaw", new.data,
        sample.id = object@sample.id,
        assembly = object@assembly, context = object@context,
        resolution = "region"
      )
    } else {
      # catch additional args
      args <- list(...)

      if (!("dbdir" %in% names(args))) {
        dbdir <- methylKit:::.check.dbdir(getwd())
      } else {
        dbdir <- methylKit:::.check.dbdir(args$dbdir)
      }
      #                         if(!( "dbtype" %in% names(args) ) ){
      #                           dbtype <- "tabix"
      #                         } else { dbtype <- args$dbtype }
      if (!("suffix" %in% names(args))) {
        suffix <- paste0("_", "regions")
      } else {
        suffix <- args$suffix
        suffix <- paste0("_", suffix)
      }

      # create methylRawDB
      obj <- methylKit:::makeMethylRawDB(
        df = new.data, dbpath = dbdir, dbtype = "tabix", sample.id = paste0(object@sample.id, suffix),
        assembly = object@assembly, context = object@context, resolution = "region"
      )
      obj@sample.id <- object@sample.id

      obj
    }
  }
)

#############
# functions #
#############

# pad regions s.t. they are 2*bp_to_pad long
pad_ranges_from_center <- function(granges_obj, bp_to_pad) {
  mid <- floor((start(granges_obj) + end(granges_obj)) / 2)
  mid_reigion_pad <- IRanges(start = (mid - bp_to_pad), end = (mid + bp_to_pad - 1))
  # -1 because it seems gr are 1 based for start and 0 based for end

  return(GRanges(seqnames = seqnames(granges_obj), ranges = mid_reigion_pad))
}

#' @description  read methylation calls/coverage for a single sample using the methylKit package
#' @param file_path path to the methylation calls file
#' @param pipeline pipeline used to generate methylation calls (default "bismarkCoverage")
#' @param sample_id sample ID as a character string
#' @param treatment treatment condition as a numeric value
#' @param assembly genome assembly used to generate methylation calls (default "mm10")
#' @param mincov minimum coverage required to include a CpG site (default 1)
#' @param context context of cytosines to include (default "CpG")
read_bismark_cov <- function(file_path,
                             sample_id,
                             assembly = "mm10",
                             mincov = 1,
                             context = "CpG",
                             pipeline = "bismarkCoverage") {
  methyl_raw <- methRead(file_path,
    sample.id = sample_id,
    assembly = assembly,
    pipeline = pipeline,
    header = TRUE,
    context = context,
    mincov = mincov
  )
  return(methyl_raw)
}

display_progress_bar <- function(current, total, bar_length = 50) {
  percent_done <- (current / total) * 100
  num_equals <- round((current / total) * bar_length)
  bar <- paste0(
    paste0(rep("=", num_equals), collapse = ""),
    paste0(rep(" ", bar_length - num_equals), collapse = "")
  )
  cat(sprintf("\r[%s] %d%% (%d/%d regions)", bar, round(percent_done), current, total))
  flush(stdout()) # Ensure the output is flushed to the console
}

process_single_region <- function(region_tiles, sample, g.meth, region_min_cov, thresh) {
  # 1. regionCounts
  # Note: regionCounts returns a methylRaw object
  counts <- regionCounts(sample, region_tiles, cov.bases = 0, strand.aware = FALSE, g.meth = g.meth)
  
  # 2. Filter low coverage
  # Direct manipulation for speed
  data <- getData(counts)
  
  # If no data (shouldn't happen with cov.bases=0 and our regionCounts mod), return NULL
  if (nrow(data) == 0) return(NULL)
  
  # Set numCs to NA where coverage is low
  data$numCs[data$coverage < region_min_cov] <- NA
  
  # 3. Check validity
  site_data <- data$numCs
  # Calculate the proportion of covered windows
  covered_proportion <- 1 - mean(is.na(site_data))
  is_valid <- covered_proportion >= thresh
  
  # 4. Calculate methylation
  # meth = 100 * numCs / (numCs + numTs)
  # This will be NA where numCs is NA
  meth_values <- 100 * data$numCs / (data$numCs + data$numTs)
  
  # 5. Prepare return object
  res <- list(
    meth_values = meth_values,
    is_valid = is_valid
  )
  
  if (is_valid) {
    # Reconstruct the region GRanges
    res$region_gr <- GRanges(
      seqnames = as.character(data$chr[1]),
      ranges = IRanges(
        start = min(data$start),
        end   = max(data$end)
      ),
      strand = as.character(data$strand[1])
    )
  }
  
  return(res)
}


parse_cli_args <- function() {
  p <- arg_parser("Average methylation by position in genomic regions")
  p <- add_argument(p, "--regions_list_file", help = "path to csv file with header including one set of regions per line: regions_name,regions_path")
  p <- add_argument(p, "--regions", help = "path to regions file - can be bed file or GRanges object saved as rds")
  p <- add_argument(p, "--regions_name", help = "name of the given regions group. used for logging and for naming the item in the results list")
  p <- add_argument(p, "--results_name", help = "results are saved as ${results_name}.rds")
  p <- add_argument(p, "--sample_name", help = "name of the sample we are calculating the methylation for")
  p <- add_argument(p, "--sample_path", help = "path to sample file - can be methylRawDB bgz or bismark coverage file")
  p <- add_argument(p, "--output_dir", help = "output directory", default = ".")
  p <- add_argument(p, "--bp_to_pad", help = "number of base pairs to pad the regions", type = "integer", default = 2500)
  p <- add_argument(p, "--tile_width", help = "width of the sliding windows ", type = "integer", default = 500)
  p <- add_argument(p, "--tile_step", help = "step of the sliding windows", type = "integer", default = 10)
  p <- add_argument(p, "--min_cov", help = "minimum coverage for a tile to be included in the average", type = "integer", default = 5L)
  p <- add_argument(p, "--mc_cores", help = "number of cores to use for parallel processing (test on one sample used ~3GB per core)", type = "integer", default = 1L)
  p <- add_argument(p, "--chunk_size", help = "number of regions to process in one batch to control memory usage", type = "integer", default = 500L)
  p <- add_argument(p, "--assembly", help = "genome assembly", default = "hg38")
  p <- add_argument(p, "--min_covered_windows", help = "min fraction of covered windows in a site to be considered valid 0.0-1.0", type = "double", default = 0.5)
  p <- add_argument(p, "--valid_sites_bed_name", help = "saves the valid sites to this file", )

  args <- parse_args(p)
  if (is.na(args$regions) && is.na(args$regions_list_file)) {
    stop("Either --regions or --regions_list_file must be provided.")
  }
  if ((!is.na(args$regions)) && is.na(args$regions_name)) {
    stop("--regions_name must be provided when using --regions")
  }
  if (is.na(args$sample_path)) {
    stop("--sample_path must be provided.")
  }
  if (is.na(args$sample_name)) {
    stop("--sample_name must be provided.")
  }
  if (is.na(args$results_name)) {
    stop("--results_name must be provided, results are saved as ${results_name}.rds")
  }
  return(args)
}

process_cli_args <- function(args) {
  is_rds <- function(file_path) {
    return(endsWith(file_path, ".rds"))
  }
  is_bed <- function(file_path) {
    return(endsWith(file_path, ".bed"))
  }
  is_methyl_db <- function(file_path) {
    return(endsWith(file_path, ".bgz"))
  }
  is_bismark <- function(file_path) {
    return(endsWith(file_path, ".cov.gz") || endsWith(file_path, ".bedGraph")) #include bedGraph from methyldackel
    #note that for methyldakel there is a header, not for bismark. for now I set the header=TRUE
    #TODO: if methyldackel is used, add a parameter to the function to set header=TRUE/FALSE
  }

  read_regions_file <- function(file_path) {
    if (is_rds(file_path)) {
      regions <- readRDS(file_path)
    } else if (is_bed(file_path)) {
      regions <- import(file_path)
    } else {
      stop("Regions file must be a bed file or a GRanges object saved as rds")
    }
    return(regions)
  }

  read_meth_file <- function(file_path, sample_id, assembly = "hg19") {
    if (is_methyl_db(file_path)) {
      # load sample to memory as a methylRaw object
      sample <- readMethylDB(file_path)
      sample <- sample[]
    } else if (is_bismark(file_path)) {
      sample <- read_bismark_cov(file_path, sample_id, assembly)
    } else {
      stop("Sample file must be a methylRawDB bgz or a bismark coverage file")
    }
    return(sample)
  }

  # read regions or regions_list_file
  if (is.na(args$regions)) {
    regions_list <- read.csv(args$regions_list_file, header = TRUE)
    cat("regions_list:\n") #debug
    print(regions_list) #debug
    # create a list of regions
    regions_names <- regions_list$regions_name
    regions_list <- lapply(regions_list$regions_path, function(regions_path) {
      cat(sprintf("reading regions file: %s\n", regions_path)) #debug
      read_regions_file(regions_path)
    })
    names(regions_list) <- regions_names

#debug
  # regions_list <- read.csv("/home/bengst/OneDriveTechnion/ido_goldstein_fasting/avg_meth_in_regions/regions_list.csv", header = TRUE)


  } else {
    regions_list <- list(read_regions_file(args$regions))
    names(regions_list) <- as.character(args$regions_name)
    print(regions_list) #debug
  }

  # debug: print a line to see if the function read_meth_file is called correctly
  cat(sprintf("debug: read_meth_file(%s, %s, %s)\n", args$sample_path, args$sample_name, args$assembly)) # debug
  sample <- read_meth_file(args$sample_path, sample_id = args$sample_name, assembly = args$assembly)
  print_mem_usage("After loading sample")
  cat("debug: class(sample) = ", class(sample), "\n", sep = "") # debug

  return(list(
    "regions_list" = regions_list,
    "sample" = sample,
    "bp_to_pad" = args$bp_to_pad,
    "tile_width" = args$tile_width,
    "tile_step" = args$tile_step,
    "min_cov" = args$min_cov,
    "mc_cores" = args$mc_cores,
    "chunk_size" = args$chunk_size
  ))
}

print_args <- function(args) {
  cat("Arguments:\n")
  if (!is.na(args$regions_list)) {
    cat("regions_list_file: ", args$regions_list_file, "\n")
    cat("regions_list:\n")
    for (i in seq_along(args$regions_list)) {
      cat(sprintf("  %s: %s\n", names(args$regions_list)[i], args$regions_list[[i]]))
    }
  } else {
    cat(sprintf("regions: %s\n", args$regions))
  }
  cat(sprintf("sample_path: %s\n", args$sample_path))
  cat(sprintf("output_dir: %s\n", args$output_dir))
  cat(sprintf("name: %s\n", args$name))
  cat(sprintf("bp_to_pad: %d\n", args$bp_to_pad))
  cat(sprintf("tile_width: %d\n", args$tile_width))
  cat(sprintf("tile_step: %d\n", args$tile_step))
  cat(sprintf("min_cov: %d\n", args$min_cov))
  cat(sprintf("mc_cores: %d\n", args$mc_cores))
  cat(sprintf("chunk_size: %d\n", args$chunk_size))
  cat(sprintf("assembly: %s\n", args$assembly))
  cat(sprintf("min_covered_windows: %f\n", args$min_covered_windows))
  cat(sprintf("valid_sites_bed_name: %s\n", args$valid_sites_bed_name))
}


main <- function(regions_list, sample, bp_to_pad, tile_width, tile_step,
                 region_min_cov, mc_cores = 1, thresh, out, valid_sites_bed_name, chunk_size = 500) {
  
  # Pre-convert sample to GRanges once
  cat("Converting sample to GRanges...\n")
  g.meth <- as(sample, "GRanges")
  strand(g.meth) <- "*" 
  print_mem_usage("After converting sample to GRanges")

  result <- lapply(names(regions_list), function(regions_name) {
    print_mem_usage(sprintf("Start processing %s", regions_name))
    regions <- regions_list[[regions_name]]
    num_regions <- length(regions)
    # print info
    cat(sprintf("Processing %s\n", regions_name))
    cat(sprintf("Number of regions: %d\n", num_regions))
    cat(sprintf(
      "padding and tiling regions:\n bp_to_pad: %d\n tile_width: %d\n tile_step: %d\n",
      bp_to_pad, tile_width, tile_step
    ))

    # pad and make sliding window tiles
    regions_pad <- pad_ranges_from_center(regions, bp_to_pad)
    regions_tiled <- slidingWindows(regions_pad, tile_width, tile_step)
    print_mem_usage("After slidingWindows")

    # Initialize accumulators
    if (length(regions_tiled) > 0) {
        num_tiles <- length(regions_tiled[[1]])
    } else {
        num_tiles <- 0
    }
    
    sum_methylation_vector <- numeric(num_tiles)
    valid_observation_count_vector <- numeric(num_tiles)
    valid_regions_gr_list <- list() 

    cat("Processing regions in chunks...\n")
    
    # Chunked execution
    num_chunks <- ceiling(num_regions / chunk_size)
    
    for (i in seq_len(num_chunks)) {
        start_idx <- (i - 1) * chunk_size + 1
        end_idx <- min(i * chunk_size, num_regions)
        
        cat(sprintf("  Chunk %d/%d (regions %d-%d)\n", i, num_chunks, start_idx, end_idx))
        
        chunk_indices <- start_idx:end_idx
        chunk_regions <- regions_tiled[chunk_indices]
        
        # Process chunk in parallel
        chunk_results <- bplapply(chunk_regions, function(region_tiles) {
            process_single_region(region_tiles, sample, g.meth, region_min_cov, thresh)
        }, BPPARAM = MulticoreParam(workers = mc_cores))
        
        # Aggregate results
        for (res in chunk_results) {
            if (!is.null(res)) {
                # Add methylation values (treat NA as 0 for sum, but don't count them)
                meth_vals <- res$meth_values
                
                # Identify non-NA values
                not_na_idx <- !is.na(meth_vals)
                
                # Update sum (replace NA with 0 for addition)
                meth_vals_clean <- meth_vals
                meth_vals_clean[is.na(meth_vals_clean)] <- 0
                sum_methylation_vector <- sum_methylation_vector + meth_vals_clean
                
                # Update counts
                valid_observation_count_vector <- valid_observation_count_vector + as.numeric(not_na_idx)
                
                # Store valid region if requested and valid
                if (!is.na(valid_sites_bed_name) && res$is_valid) {
                    valid_regions_gr_list[[length(valid_regions_gr_list) + 1]] <- res$region_gr
                }
            }
        }
        
        # Cleanup
        rm(chunk_results)
        gc()
    }
    
    print_mem_usage("After processing all chunks")

    # Export valid sites if requested
    if (!is.na(valid_sites_bed_name) && length(valid_regions_gr_list) > 0) {
        cat("Exporting valid sites...\n")
        valid_regions_gr <- do.call(c, valid_regions_gr_list)
        export.bed(valid_regions_gr, valid_sites_bed_name)
    }

    # Calculate average
    cat("Calculating final average...\n")
    avg_meth_by_pos <- sum_methylation_vector / valid_observation_count_vector
    
    dist_from_center <- seq(-bp_to_pad, (bp_to_pad - 1), length.out = length(avg_meth_by_pos))

    return(list("avg_meth_by_pos" = avg_meth_by_pos, "dist_from_center" = dist_from_center))
  })
  names(result) <- names(regions_list) #needs testing (18.10.24)
  return(result)
}


########################
# run from command line#
########################
args <- parse_cli_args()
print_args(args)
args_main <- process_cli_args(args)
result <- main(
  regions_list = args_main$regions_list,
  sample = args_main$sample,
  bp_to_pad = args_main$bp_to_pad,
  tile_width = args_main$tile_width,
  tile_step = args_main$tile_step,
  region_min_cov = args_main$min_cov,
  mc_cores = args_main$mc_cores,
  thresh = args$min_covered_windows,
  out = args$output_dir,
  valid_sites_bed_name = args$valid_sites_bed_name,
  chunk_size = args_main$chunk_size
)
# if output dir does not exist, create it
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir)
}
cat(sprintf("saving file: %s/%s.rds", args$output_dir, args$results_name))
saveRDS(result, file = sprintf("%s/%s.rds", args$output_dir, args$results_name))
quit(status = 0)



############
# comments #
############

# [1]
# order of regions may change with parallelization, no problem since we take the average.
# order of tiles in each region is preserved (a region's tiles set is handeled by one task).

# TODO: look at method for regions = GRangesList and see if they are doing it faster than lapply,
#  if so take that and change signature to accept (comressed/simple)GRangesList - because I can't seem
#  to be able to produce a GRangesList from a list of GRanges, maybe deprecated?

# TODO: test how efficient the parallelization is by timing the script using 1,10 cores.

# TODO: may be faster to coerce list of GRanges to GRangesList and run regionCounts once
# does not work - I get compresssedGRangesList or simpleGRangesList instead of GRangesList
# #convert list of GRanges to GRangesList or set method for compressedGRangesList
# regions_tiled_grlist <- GRangesList(regions_tiled, compress = FALSE)
#
# time_regioncounts_grlist <- system.time({
#   counts2 <- regionCounts(sample, as(regions_tiled, "GRangesList"))
# })



# [2]
# could have done this in regionCounts method above but
# decided to leave it as close to the original as possible
# also I think the logic flow is less confusing this way


# [4]
# Tested 3 methods for replacing numCs with NA in low coverage tiles
# fastest was lapply with direct manipulation of the data.
# below are results and then the functions
#
# results
# Processing N-ARBS
# Number of regions: 50
# padding and tiling regions:
#  bp_to_pad: 2500
#  tile_width: 500
#  tile_step: 10
# Counting number of Cs and Ts in each tile
# Replacing numCs with NA in tiles with coverage < 5. Calculated methylation will be NA
# loop: 12.649000
# lapply_1: 0.014000
# lapply_2: 0.020000
# Method 1 and Method 2 give the same result
# Method 1 and Method 3 give the same result
# Method 2 and Method 3 give the same result
# Calculating average methylation per position
#
# Processing T-ARBS
# Number of regions: 50
# padding and tiling regions:
#  bp_to_pad: 2500
#  tile_width: 500
#  tile_step: 10
# Counting number of Cs and Ts in each tile
# Replacing numCs with NA in tiles with coverage < 5. Calculated methylation will be NA
# loop: 13.206000
# lapply_1: 0.011000
# lapply_2: 0.026000
# Method 1 and Method 2 give the same result
# Method 1 and Method 3 give the same result
# Method 2 and Method 3 give the same result
# Calculating average methylation per position
#
#
#
# functions
# using looop to replace numCs with NA in low coverage tiles
# numcs_na_if_low_coverage_loop <- function(counts, region_min_cov) {
#   for (i in seq_along(counts)) {
#     for (j in seq_len(nrow(counts[[i]]))) {
#       if (counts[[i]][j]$coverage < region_min_cov) {
#         counts[[i]][j, "numCs"] <- NA
#       }
#     }
#   }
#   return(counts)
# }
#
# # using laapply to replace numCs with NA in low coverage tiles
# numcs_na_if_low_coverage_lapply_1 <- function(counts, region_min_cov) {
#   counts <- lapply(counts, function(x) {
#     x[x$coverage < region_min_cov, "numCs"] <- NA
#     return(x)
#   })
#   return(counts)
# }
#
# numcs_na_if_low_coverage_lapply_2 <- function(counts, region_min_cov) {
#   counts <- lapply(counts, function(methylRawObj) {
#     data <- getData(methylRawObj)
#     data$numCs[data$coverage < region_min_cov] <- NA
#     return(new("methylRaw", data,
#       sample.id = methylRawObj@sample.id,
#       assembly = methylRawObj@assembly,
#       context = methylRawObj@context,
#       resolut########################
# run from command line#
########################
args <- parse_cli_args()
print_args(args)
args_mion = "region"
#     ))
#   })
#   return(counts)
# }
# ./14_clinical_samp_skvortsova2019/methylseq/mbias_trim/work
#
# testing
# measure and print run time of different functions for this task
# time_numcs_na_if_low_coverage_loop <- system.time({
#   counts1 <- numcs_na_if_low_coverage_loop(counts, region_min_cov)
# })
#
# time_numcs_na_if_low_coverage_lapply_1 <- system.time({
#   counts2 <- numcs_na_if_low_coverage_lapply_1(counts, region_min_cov)
# })
#
# time_numcs_na_if_low_coverage_lapply_2 <- system.time({
#   counts3 <- numcs_na_if_low_coverage_lapply_2(counts, region_min_cov)
# })
#
# cat(sprintf("loop: %f\n", time_numcs_na_if_low_coverage_loop["elapsed"]))
# cat(sprintf("lapply_1: %f\n", time_numcs_na_if_low_coverage_lapply_1["elapsed"]))
# cat(sprintf("lapply_2: %f\n", time_numcs_na_if_low_coverage_lapply_2["elapsed"]))
#
# if (identical(counts1, counts2)) {
#   cat("Method 1 and Method 2 give the same result\n")
# } else {
#   cat("Method 1 and Method 2 give different results\n")
# }
#
# if (identical(counts1, counts3)) {
#   cat("Method 1 and Method 3 give the same result\n")
# } else {
#   cat("Method 1 and Method 3 give different results\n")
# }
#
# if (identical(counts2, counts3)) {
#   cat("Method 2 and Method 3 give the same result\n")
# } else {
#   cat("Method 2 and Method 3 give different results\n")
# }
