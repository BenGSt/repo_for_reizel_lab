#!/usr/bin/env Rscript
library(argparser)
library(stringr)
library(magrittr)
library("pheatmap")


main = function(scores_bed_file, sample_names=NULL, include_samples_by_name=NULL,
                title="Methylation Scores", output_file="heatmap.png")
{
  meth_scores = read.delim(scores_bed_file)
  if (! is.null(sample_names))
    colnames(meth_scores)[4:ncol(meth_scores)] = sample_names 
  
  if (! is.null(include_samples_by_name))
  {
    include_columns = colnames(meth_scores)  %in% include_samples_by_name 
    meth_scores = cbind(meth_scores[,1:3], meth_scores[,include_columns])
  }

  meth_matrix = meth_scores[4:ncol(meth_scores)] %>% as.matrix()
  
  
  col <- colorRampPalette(c("blue", "yellow"))(256)
  png(file=output_file, width=1000, height=1000, res=150)
  pheatmap(meth_matrix, cutree_cols=2, color=col, cluster_rows=FALSE, main=title)
  dev.off()
}


# Create a parser
p <- arg_parser("Make DMRs heatmap")

# Add command line arguments
p <- add_argument(p, "--scores_bed_file", short="-m", help="bed file with methylation scores for the samples")
p <- add_argument(p, "--sample_names", short="-s", help="string with the names of the samples seperated by \"-\" (optional)", default=NULL)
p <- add_argument(p, "--include_samples_by_name", short="-i", help="names of the samples to include in the heatmap (subset of sample_names)", default=NULL)
p <- add_argument(p, "--output_file", short="-o", help="path to output png", default = "heatmap.png")
p <- add_argument(p, "--title", short="-t", help="plot title", default = "Methylation Scores Heatmap")

# Parse the command line arguments
argv <- parse_args(p)


if (! is.na(argv$sample_names)){
  sample_names = strsplit(argv$sample_names,'-')[[1]]
} else {
  sample_names = NULL
  }

if (! is.na(argv$include_samples_by_name)) {
  include_samples_by_name = strsplit(argv$include_samples_by_name,'-')[[1]]
} else {
  include_samples_by_name = NULL
}

main(argv$scores_bed_file, sample_names, include_samples_by_name, str_replace(argv$title,"\\\\n","\n"), argv$output_file)
