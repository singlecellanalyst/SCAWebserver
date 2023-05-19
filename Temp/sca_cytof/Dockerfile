FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scacytof/
COPY ui.R /srv/shiny-server/scacytof/
COPY DB/ /srv/shiny-server/scacytof/DB 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y openjdk-17-jre-headless && \
    apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(cran <- c("shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","remotes","ggplot2","reshape2","gplots","plot3D","Rtsne","wesanderson","viridis","webshot","plotly","lattice","latticeExtra","ggpubr","plyr","akmedoids","umap","patchwork","dplyr","ggridges","ggthemes","gridExtra","ggrepel","scales","xlsx","BiocManager","devtools","DT","Seurat","rJava"))'

RUN Rscript -e 'BiocManager::install(c("flowCore","limma","Biobase","FlowSOM","ConsensusClusterPlus","ComplexHeatmap"))'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scacytof

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
