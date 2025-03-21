FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scasccnv/
     COPY ui.R /srv/shiny-server/scasccnv/
     COPY DB/ /srv/shiny-server/scasccnv/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y default-jdk && \
apt-get install -y r-cran-rjava && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk libglu1 && \
rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("shiny","BiocManager","devtools","remotes","dplyr","ggplot2","rgl","R.utils"))'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(cran <- c("Seurat","patchwork","ggridges","ggthemes","gridExtra","reshape2","ggrepel","scales","xlsx","shinythemes","shinyWidgets","shinyscreenshot","shinyjs","pheatmap","cowplot","ape","data.table","intervalaverage","RColorBrewer","adegenet","fastcluster","umap"))'

RUN Rscript -e 'BiocManager::install(c("GenomicRanges","ComplexHeatmap","ggtree"),ask = F)'

RUN Rscript -e 'devtools::install_github(repo = "jmw86069/jamba",force = T)'
RUN Rscript -e 'devtools::install_github(repo = "jmw86069/splicejam",force = T)'

RUN Rscript -e 'packages <- c("shiny","BiocManager","devtools","remotes","dplyr","ggplot2","rgl","R.utils","Seurat","patchwork","ggridges","ggthemes","gridExtra","reshape2","ggrepel","scales","xlsx","shinythemes","shinyWidgets","shinyscreenshot","shinyjs","pheatmap","cowplot","ape","data.table","intervalaverage","RColorBrewer","adegenet","GenomicRanges","ComplexHeatmap","ggtree","jamba","splicejam","fastcluster"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/scasccnv

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]