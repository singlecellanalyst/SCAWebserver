FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scascimmune/
  COPY ui.R /srv/shiny-server/scascimmune/
  COPY DB/ /srv/shiny-server/scascimmune/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y default-jdk && \
apt-get install -y r-cran-rjava && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk && \
rm -rf /var/lib/apt/lists/*
  
RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("igraph","immunarch","shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","webshot","ggthemes","plot3D","ggalluvial","dplyr","alluvial","harmony","ggrepel","RColorBrewer","circlize","scales","gridExtra","hdf5r","ggraph","assertthat","DT","ggplot2","patchwork","ggridges","reshape2","xlsx","BiocManager","httr","rJava","Seurat","devtools","remotes"))'

RUN Rscript -e 'BiocManager::install(c("monocle","SingleCellExperiment","powerTCR","SingleR","ComplexHeatmap"), ask = F)'
RUN Rscript -e 'devtools::install_github("ncborcherding/scRepertoire", force = T)'

RUN Rscript -e 'packages <- c("igraph","immunarch","shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","webshot","ggthemes","plot3D","ggalluvial","dplyr","alluvial","harmony","ggrepel","RColorBrewer","circlize","scales","gridExtra","hdf5r","ggraph","assertthat","DT","ggplot2","patchwork","ggridges","reshape2","xlsx","BiocManager","httr","rJava","Seurat","devtools","remotes","monocle","SingleCellExperiment","powerTCR","SingleR","ComplexHeatmap","scRepertoire"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/scascimmune

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]