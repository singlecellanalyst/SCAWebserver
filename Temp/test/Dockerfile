FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/test/
  COPY ui.R /srv/shiny-server/test/
  COPY DB/ /srv/shiny-server/test/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y openjdk-17-jre-headless && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk && \
rm -rf /var/lib/apt/lists/*

RUN which java
ENV LD_LIBRARY_PATH=/usr/lib/jvm/java-17-openjdk-amd64/

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("igraph","shiny","shinyjs","BiocManager","httr","rJava","xlsx","Seurat","devtools","remotes"))'

RUN Rscript -e 'packages <- c("igraph","shiny","shinyjs","BiocManager","httr","rJava","xlsx","Seurat","devtools","remotes"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/test

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]