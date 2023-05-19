FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scaflow/
COPY ui.R /srv/shiny-server/scaflow/
COPY DB/ /srv/shiny-server/scaflow/DB 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y openjdk-17-jre-headless && \
    apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(cran <- c("xlsx","shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","remotes","robustbase","SDMTools","rhandsontable","gridExtra","reshape2","ggridges","patchwork","gplots","plot3D","DT","umap","viridis","cowplot","akmedoids","gridGraphics","ggplot2","dplyr","ggthemes","ggrepel","scales","BiocManager","devtools","remotes","rJava")'

RUN Rscript -e 'BiocManager::install(c("Seurat","cytolib","flowCore","ggcyto","flowStats","openCyto","flowViz","limma","Biobase","cytofWorkflow","FlowSOM","ConsensusClusterPlus","flowWorkspace","flowDensity"))'

RUN Rscript -e 'devtools::install_github("ParkerICI/vite")'
RUN Rscript -e 'devtools::install_github(repo = "zunderlab/FLOWMAP")'
RUN Rscript -e 'devtools::install_github("DillonHammill/CytoExploreRData")'
RUN Rscript -e 'devtools::install_github("DillonHammill/CytoExploreR")'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scaflow

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
