FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scascatac/
COPY ui.R /srv/shiny-server/scascatac/
COPY DB/ /srv/shiny-server/scascatac/DB 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y openjdk-17-jre-headless && \
    apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(cran <- c("remotes","devtools","BiocManager","R.utils"))'
RUN Rscript -e 'devtools::install_github("thomasp85/patchwork")")'
RUN Rscript -e 'BiocManager::install(c("AnnotationFilter","BiocGenerics","GenomicFeatures","GenomicRanges","IRanges","Rsamtools","S4Vectors","ggbio","motifmatchr","AnnotationDbi","Seurat","GenomeInfoDb","JASPAR2020","TFBSTools","EnsDb.Hsapiens.v86","BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19"))'

RUN Rscript -e 'install.packages("Signac")'

RUN Rscript -e 'install.packages(cran <- c("shiny","shinythemes","ggplot2","DT","shinyWidgets","shinyscreenshot","patchwork","dplyr","ggridges","ggthemes","gridExtra","reshape2","ggrepel","scales","xlsx","shinyjs","shinyFiles","shinydashboard","shinyalert","webshot","gpubr","cowplot","hdf5r","qlcMatrix","ggseqlogo","ComplexHeatmap","rJava"))'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scascatac

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
