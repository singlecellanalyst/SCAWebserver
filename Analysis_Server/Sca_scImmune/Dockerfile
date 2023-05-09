FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libmagick++-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","webshot","ggthemes","plot3D","ggalluvial","dplyr","alluvial","harmony","ggrepel","RColorBrewer","immunarch","circlize","scales","gridExtra","hdf5r","ggraph","assertthat","DT","ggplot2","patchwork","ggridges","reshape2","xlsx"))'
RUN Rscript -e 'BiocManager::install(c("Seurat","monocle","SingleCellExperiment","scRepertoire","powerTCR","SingleR","ComplexHeatmap")'

COPY /DB/ /srv/shiny-server/app/DB 
COPY global.R /srv/shiny-server/app/
COPY server.R /srv/shiny-server/app/
COPY ui.R /srv/shiny-server/app/

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/app

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
