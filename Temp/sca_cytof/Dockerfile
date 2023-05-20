FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scacytof/
COPY ui.R /srv/shiny-server/scacytof/
COPY DB/ /srv/shiny-server/scacytof/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y default-jdk && \
apt-get install -y r-cran-rjava && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk libglu1 && \
rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("shiny","BiocManager","devtools","remotes","kml","dplyr","signal","Hmisc","ggplot2","rgl"))'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages("https://cran.r-project.org/src/contrib/Archive/clusterCrit/clusterCrit_1.2.8.tar.gz")'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages("https://cran.r-project.org/src/contrib/Archive/akmedoids/akmedoids_1.3.0.tar.gz")'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(cran <- c("shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","reshape2","gplots","plot3D","Rtsne","wesanderson","viridis","webshot","plotly","lattice","latticeExtra","ggpubr","plyr","umap","patchwork","ggridges","ggthemes","gridExtra","ggrepel","scales","xlsx","DT","Seurat","rJava"))'

RUN Rscript -e 'BiocManager::install(c("flowCore","limma","Biobase","FlowSOM","ConsensusClusterPlus","ComplexHeatmap"), ask = F)'
RUN Rscript -e 'packages <- c("shiny","shinyjs","shinyWidgets","shinyscreenshot","shinythemes","shinyFiles","shinydashboard","shinyalert","remotes","ggplot2","reshape2","gplots","plot3D","Rtsne","wesanderson","viridis","webshot","plotly","lattice","latticeExtra","ggpubr","plyr","umap","patchwork","dplyr","ggridges","ggthemes","gridExtra","ggrepel","scales","xlsx","devtools","DT","Seurat","rJava","flowCore","limma","Biobase","FlowSOM","ConsensusClusterPlus","ComplexHeatmap"); lapply(packages, require, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scacytof

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
