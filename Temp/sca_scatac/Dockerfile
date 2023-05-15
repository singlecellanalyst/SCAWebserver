FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/app/
COPY ui.R /srv/shiny-server/app/
COPY DB/ /srv/shiny-server/app/DB 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(cran <- c("shiny","shinythemes","ggplot2","DT","shinyWidgets","shinyscreenshot","patchwork","dplyr","ggridges","ggthemes","gridExtra","reshape2","ggrepel","scales","xlsx","BiocManager","shinyjs","shinyFiles","shinydashboard","shinyalert","webshot","Signac","ggpubr","cowplot","hdf5r","qlcMatrix","ggseqlogo","remotes","devtools"))'

RUN Rscript -e 'BiocManager::install(c("Seurat","GenomeInfoDb","motifmatchr","JASPAR2020","TFBSTools","EnsDb.Hsapiens.v86","BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19"))'

RUN Rscript -e 'remotes::install_github('satijalab/seurat-wrappers')'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/app

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
