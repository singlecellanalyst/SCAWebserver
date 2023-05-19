FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/test/
  COPY ui.R /srv/shiny-server/test/
  COPY DB/ /srv/shiny-server/test/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev && \
rm -rf /var/lib/apt/lists/*
  
RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("devtools","plot3D","igraph","shiny","BiocManager"))'
RUN Rscript -e 'packages <- c("plot3D","igraph","shiny","scRepertoire"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/test

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]