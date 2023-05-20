FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/test/
  COPY ui.R /srv/shiny-server/test/
  COPY DB/ /srv/shiny-server/test/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y default-jdk && \
apt-get install -y r-cran-rjava && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk libglu1 && \
rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("shiny","BiocManager","devtools","remotes","ggplot2"))'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages("https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz")'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(cran <- c("shiny","shinyjs","robustbase","rhandsontable","gridExtra","reshape2","ggridges","patchwork","gplots","plot3D","umap","viridis","cowplot","gridGraphics","ggthemes","ggrepel","scales","rJava"))'

RUN Rscript -e 'BiocManager::install(c("cytolib","flowCore","ggcyto","flowStats","openCyto","flowViz","limma","Biobase","cytofWorkflow","FlowSOM","ConsensusClusterPlus","flowWorkspace","flowDensity"),ask = F)'

RUN Rscript -e 'devtools::install_github("ParkerICI/vite",force = T)'
RUN Rscript -e 'devtools::install_github(repo = "zunderlab/FLOWMAP",force = T)'
RUN Rscript -e 'devtools::install_github("DillonHammill/CytoExploreRData",force = T)'
RUN Rscript -e 'devtools::install_github("DillonHammill/CytoExploreR",force = T)'

RUN Rscript -e 'packages <- c("akmedoids","shiny"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/test

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]