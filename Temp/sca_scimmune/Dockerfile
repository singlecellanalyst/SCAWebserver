FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scascimmune/
COPY ui.R /srv/shiny-server/scascimmune/
COPY DB/ /srv/shiny-server/scascimmune/DB 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("igraph","immunarch","shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","webshot","ggthemes","plot3D","ggalluvial","dplyr","alluvial","harmony","ggrepel","RColorBrewer","circlize","scales","gridExtra","hdf5r","ggraph","assertthat","DT","ggplot2","patchwork","ggridges","reshape2","xlsx","BiocManager"))'

RUN Rscript -e 'BiocManager::install(c("Seurat","monocle","SingleCellExperiment","scRepertoire","powerTCR","SingleR","ComplexHeatmap"))'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scascimmune

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]

