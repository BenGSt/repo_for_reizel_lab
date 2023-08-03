#!/usr/bin/env Rscript
suppressMessages(library(DESeq2))
suppressMessages(library(ReportingTools))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library(argparser))

GTF_USED_BY_HTSEQ = 'ENSEMBL' # keytype argument for "org.Hs.eg.db"
N_GENES_REPORT = 500 # number of most significantly differentialy expressed genes to report in the html report

#' MakeSampleTable
#' Make a table summarizing all the htseq-count output files
#' @param htseq_out_dir directory containing all the htseq-count output files of
#'                    the samples in the experiment
#'
#' @return sample_table (dataframe)
MakeSampleTable = function(htseq_out_dir)
{
  sample_files <- list.files(htseq_out_dir)
  sample_conditions <- sub("(^.*)[0-9].*$", "\\1", sample_files)
  sample_names <- sub("(^.*[0-9]).*$", "\\1", sample_files)
  sample_table <- data.frame(sample_name = sample_names,
                             file = sample_files,
                             condition = sample_conditions)

  # change condition data type to "factor" object
  sample_table$condition <- factor(sample_table$condition)

  return(sample_table)
}


#' buildDESeqDataSet
#'
#' @param sample_table dataframe returned from makeSampleTable()
#' @param htseq_out_dir directory containing all the htseq-count output files of
#'                    the samples in the experiment
#'
#' @return DESeqDataSet object
buildDESeqDataSet = function(sample_table, htseq_out_dir)
{
  dds <- DESeqDataSetFromHTSeqCount(sample_table, htseq_out_dir, ~condition)

  # pre-filtering to keep the dds object's size down and speed things up.
  # keep only rows that have at least 10 reads total.
  # Note that more strict filtering to increase power is automatically applied via
  # independent filtering on the mean of normalized counts within the results function.
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  return(dds)
}

# dirname = function(path)
# {
#   split = strsplit(path,"/")[[1]]
#   dirname = split[[length(split)]]
#   return(dirname)
# }


buildHtmlReport = function(dds, report_dir, pvalue_cutoff = 0.01, contrast = NULL)
{
  # NOTE: pvalueCutoff here refers to adjusted pvalue in report
  #   report_dirname = dirname(report_dir)
  #   setwd("../")

  colData(dds)$conditions <- dds$condition
  test = strsplit(results(dds, contrast = contrast)@
                    elementMetadata@
                    listData$description[2], ":")[[1]][2]
  report <- HTMLReport(shortName = 'deg_analysis_with_deseq2',
                       title = paste('RNA-seq analysis of differential expression using DESeq2\n',
                                     test, 'padj_cutoff=', pvalue_cutoff),
                       reportDirectory = report_dir)

  publish(dds, report, pvalueCutoff = pvalue_cutoff,
          n = N_GENES_REPORT,
          factor = colData(dds)$conditions,
          annotation.db = "annotation_db",
          keytype = GTF_USED_BY_HTSEQ,
          contrast = contrast,
          reportDir = report_dir)

  finish(report)
  on.exit(setwd(report_dir))
}


volcano = function(res, pvalue_cutoff = 10e-10, log2_fc_cutoff = 0.5, png_path, title)
{
  symbols <- mapIds(annotation_db, keys = rownames(res),
                    column = c('SYMBOL'), keytype = GTF_USED_BY_HTSEQ)
  symbols[is.na(symbols)] <- names(symbols[is.na(symbols)])
  rownames(res) <- symbols
  cap = paste("log2 fold-change cutoff: ", log2_fc_cutoff, ";  p-value cutoff: ", pvalue_cutoff)

  png(file = png_path, width = 2000, height = 1500, res = 150)
  print(
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title = title,
                    drawConnectors = T,
                    caption = cap,
                    subtitle = "",
                    pCutoff = pvalue_cutoff,
                    FCcutoff = log2_fc_cutoff,
                    pointSize = 3.0,
                    labSize = 4.0
    ) +
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}

volcano_padj = function(res, adjp_cutoff = 10e-10, log2_fc_cutoff = 0.5, png_path, title)
{
  symbols <- mapIds(annotation_db, keys = rownames(res),
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  symbols[is.na(symbols)] <- names(symbols[is.na(symbols)])
  rownames(res) <- symbols
  cap = paste("log2 fold-change cutoff: ", log2_fc_cutoff, ";  adjusted p-value cutoff: ", adjp_cutoff)

  png(file = png_path, width = 2000, height = 1500, res = 150)
  print(
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = title,
                    drawConnectors = T,
                    caption = cap,
                    subtitle = "",
                    pCutoff = adjp_cutoff,
                    FCcutoff = log2_fc_cutoff,
                    pCutoffCol = "padj",
                    pointSize = 3.0,
                    labSize = 4.0
    ) +
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}


main = function(htseq_out_dir, report_dir, padj_cutoff = 0.01, log2_fc_cutoff = 1,
                contrast = NULL, csv_path = "deg_deseq2.csv", volcano_plot_title = NULL, png_path = "volcano_plot.png")
{


  sampleTable = MakeSampleTable(htseq_out_dir)
  dds = buildDESeqDataSet(sampleTable, htseq_out_dir)
  dds <- DESeq(dds)

  #TEST: 3.8.2023 added PCA
  vsd <- vst(dds, blind = FALSE)
  png(file = "pca.png", width = 2000, height = 1500, res = 150)
  print(plotPCA(vsd))
  dev.off()

  buildHtmlReport(dds, report_dir, padj_cutoff, contrast)

  # By default the argument alpha is set to 0.1. If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value
  res <- results(dds, contrast = contrast, alpha = padj_cutoff)
  volcano_padj(res, padj_cutoff, log2_fc_cutoff, png_path, volcano_plot_title)
  deg_list = res[abs(res$log2FoldChange) >= log2_fc_cutoff &
                   !is.na(res$padj) &
                   res$padj <= padj_cutoff,]
  symbols <- mapIds(annotation_db, keys = rownames(deg_list),
                    column = c('SYMBOL'), keytype = GTF_USED_BY_HTSEQ)
  symbols[is.na(symbols)] <- names(symbols[is.na(symbols)])
  rownames(deg_list) <- symbols
  write.csv(as.data.frame(deg_list), file = csv_path)
}


p <- arg_parser("differential expression analysis with DESeq2")
p <- add_argument(p, "--htseq_output_dir", help = "directory containing the output files of htseq-count")
p <- add_argument(p, "--report_dir", help = "directory to write DEG HTML report to, as well as csv and volacano plot")
p <- add_argument(p, "--padj_cutoff", help = "the adjusted pvalue cutoff value for the volcano plot and csv file", type = "numeric")
p <- add_argument(p, "--log2_fc_cutoff", help = "the log2 fold change cutoff value for the volcano plot and csv file", type = "numeric")
p <- add_argument(p, "--contrast", help = "specifies what comparison to extract from the object to build a results table e.g. c(\"condition\",\"F\",\"I\"). see ?DESeq2::results")
p <- add_argument(p, "--volcano_plot_title", help = "the title for the volcano plot")
p <- add_argument(p, "--csv", help = "file to write csv deg list to")
p <- add_argument(p, "--genome", help = "mm10 or hg38", type = "character")
argv <- parse_args(p)

if (is.na(argv$volcano_plot_title))
{
  contrast = eval(parse(text = argv$contrast))
  argv$volcano_plot_title = paste(contrast[2], "vs.", contrast[3])
}

if (argv$genome == "mm10") {
  annotation_db <- org.Mm.eg.db
} else if (argv$genome == "hg38") {
  annotation_db <- org.Hs.eg.db
} else {
  print(sprintf("genome %s not supported", argv$genome))
  quit(status = 1)
}

main(
  htseq_out_dir = argv$htseq_output_dir,
  report_dir = argv$report_dir,
  padj_cutoff = argv$padj_cutoff,
  log2_fc_cutoff = argv$log2_fc_cutoff,
  contrast = eval(parse(text = argv$contrast)),
  csv_path = argv$csv,
  volcano_plot_title = argv$volcano_plot_title,
  png_path = "volcano_plot.png"
)


#Note:  the error "'select()' returned 1:many mapping between keys and columns" Is because of NAs returned for some keys.
# Looks like mostly novel transcripts. This is fixed in the csv file by reverting to the original ensembl id.
# In the HTML report NA will appear. could fix this if I have the time, probably not worth it