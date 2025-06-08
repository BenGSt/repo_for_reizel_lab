# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install required R packages
RUN R -e 'install.packages(c("argparser", "dplyr", "stringr"))'

# Install required Bioconductor package
RUN R -e 'BiocManager::install("methylKit")'

