FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scaspatial/
  COPY ui.R /srv/shiny-server/scaspatial/
  COPY DB/ /srv/shiny-server/scaspatial/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y default-jdk && \
apt-get install -y r-cran-rjava && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk libglu1 && \
rm -rf /var/lib/apt/lists/*
  
RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("Seurat","ggplot2","patchwork","dplyr","ggridges","ggthemes","gridExtra","reshape2","ggrepel","scales","xlsx","shiny","shinythemes","DT","shinyWidgets","shinyscreenshot","shinyjs","ggpubr","hdf5r","BiocManager","devtools"))'

RUN Rscript -e 'BiocManager::install(c("SingleR","org.Hs.eg.db","SingleCellExperiment","limma"), ask = F)'

RUN Rscript -e 'packages <- c("Seurat","ggplot2","patchwork","dplyr","ggridges","ggthemes","gridExtra","reshape2","ggrepel","scales","xlsx","shiny","shinythemes","DT","shinyWidgets","shinyscreenshot","shinyjs","ggpubr","hdf5r","SingleR","org.Hs.eg.db","SingleCellExperiment","limma"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/scaspatial

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]