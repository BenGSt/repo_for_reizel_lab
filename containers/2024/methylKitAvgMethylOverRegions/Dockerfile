# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Install required R packages
RUN R -e 'install.packages(c( \
     "rtracklayer", \
     "data.table", \
     "dplyr", \
     "parallel", \
     "argparser", \
     "IRanges"))'

# Install required Bioconductor packages
#library(GenomeInfoDb)
#library(GenomicRanges)
#library(methylKit)
RUN R -e 'BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "methylKit"))'
