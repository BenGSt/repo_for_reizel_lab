# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_19


# Install required Bioconductor package
RUN R -e 'BiocManager::install(c(\
"annotatr",\
"org.Hs.eg.db",\
"TxDb.Hsapiens.UCSC.hg19.knownGene",\
"TxDb.Hsapiens.UCSC.hg38.knownGene",\
"org.Mm.eg.db",\
"TxDb.Mmusculus.UCSC.mm10.knownGene",\
"TxDb.Mmusculus.UCSC.mm10.knownGene"\
))'

# Install required R packages
RUN R -e 'install.packages(c("dplyr", "ggplot2", "argparser"))'