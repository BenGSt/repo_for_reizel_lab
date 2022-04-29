library(stringr)
library(magrittr)
#FAH regeneration info
experiments_barcodes = c("AAGAGG", "AAGCCT", "ACCTCA", "AGCATG", "AGTGAG",
                         "GAGTCA", "GCACTA", "GGAGAA", "GTCGTA", "GTGCTT","TGGTGA")
samp_ids=c("YoungYoung", "YoungYoung", "OldOld", "YoungYoungProlong",
           "Young","YoungYoungProlong", "OldOld", "Old", "Young", "OldOld", "Old")
names(samp_ids) = experiments_barcodes

#get data
setwd("C:/Users/bengs/Nextcloud/Tzachi_bioinformatics/Fah_regeneration/heatmaps")
bed.files = list.files(pattern = ".bed")
df_names = bed.files %>%  str_replace(".bed","")
meth_scores = lapply(bed.files, read.delim)
names(meth_scores) = df_names
meth_scores = lapply(meth_scores,
                       function(x) 
                       {
                         #replace barcodes with experiment name
                         colnames(x)[4:ncol(x)] =  
                            lapply(colnames(x)[4:ncol(x)],function(y) samp_ids[y])
                         return(x)
                       }
                     )

meth_matrixes = lapply(meth_scores, function(x) x[4:ncol(x)] %>% as.matrix())
#TODO: make row names the genomic coordinates

mat_index = 3
meth_mat = meth_matrixes[[mat_index]]

###################################
###########R base heatmap##########
###################################
#scale: a character indicating if the values should be 
# centered and scaled in either the row direction or the column direction,
# or none. Allowed values are in c(“row”, “column”, “none”). Default is “row”.
# Centering is done by subtracting the means, scaling is done by dividing
# (centered) values by their standard deviations

heatmap(meth_matrixes[[4]], scale="none",hclustfun = function(x) hclust(x, method = "complete"))

library("RColorBrewer")
col <- colorRampPalette(brewer.pal(9, "Blues"))(256)
col <- colorRampPalette(c("blue", "yellow"))(256)
#column_colors = lapply((meth_mat %>% colnames() %>% str_remove(pattern = ".[0-9]")),FUN = function(x) experiment_color_function[x] ) %>% as.character()


heatmap(meth_mat, scale = "none", col =  col, Rowv=NA,
        #RowSideColors = rep(c("blue", "pink","red","green","yellow"), each = (nrow(meth_mat) / 5)),
       # ColSideColors = column_colors
       )

##########################################
#####Enhanced heat maps: heatmap.2()######
##########################################
#install.packages("gplots")
library("gplots")
col <- colorRampPalette(c("blue", "yellow"))(256)


heatmap.2(meth_mat, scale = "none", col = col,
          dendrogram = "column",labRow = FALSE,
          trace = "none", density.info = "none")

########################################
#####Pretty heat maps: pheatmap()#######
########################################
#install.packages("pheatmap")
library("pheatmap")



plot_pheatmaps_separately = function(i, directory=getwd())
{
  mat_index = i
  meth_mat = meth_matrixes[[mat_index]]
  mat_name = names(meth_matrixes)[mat_index]
  old_wd = getwd()
  setwd(directory)
  png(file=str_c(mat_name,"_heatmap.png"),width=1000,height=1000,res=150)
  pheatmap(meth_mat, cutree_cols = 2, color = col, cluster_rows=FALSE,
           main=str_c("FAH Regeneration Methylation Scores\n",str_c ("DMRs: " ,mat_name %>% str_remove("_meth_scores_dmrs"))))
  dev.off()
  setwd(old_wd)
}


plot_pheatmaps_separately_no_prolong = function(i, directory=getwd())
{
  include_columns = meth_matrixes[[mat_index]] %>% colnames  %in% c("YoungYoungProlong","YoungYoungProlong.1") %>% not()
  mat_index = i
  meth_mat = meth_matrixes[[mat_index]][,include_columns]
  mat_name = names(meth_matrixes)[mat_index]
  old_wd = getwd()
  setwd(directory)
  png(file=str_c(mat_name,"_noprolong_heatmap.png"),width=1000,height=1000,res=150)
  pheatmap(meth_mat, cutree_cols = 2, color = col, cluster_rows=FALSE,
           main=str_c("FAH Regeneration Methylation Scores\n(Exluding Prolonged Samples)\n",str_c ("DMRs: " ,mat_name %>% str_remove("_meth_scores_dmrs"))))
  dev.off()
  setwd(old_wd)
}


dir.create("plot_pheatmaps_separately")
lapply(seq(1,4), function(i) plot_pheatmaps_separately(i,directory="plot_pheatmaps_separately"))
lapply(seq(1,4), function(i) plot_pheatmaps_separately_no_prolong(i,directory="plot_pheatmaps_separately"))
# seq out of bounds intentionally - last plot comes out screwed up if not.


#######################################
####Interactive heat maps: d3heatmap() - not available in R4
######################################
# install.packages("d3heatmap")
library("d3heatmap")
d3heatmap(meth_mat, colors = "RdYlBu",
          k_col = 2 # Number of groups in columns
)

# see more at https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
