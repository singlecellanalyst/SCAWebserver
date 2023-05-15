FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scaflow/
COPY ui.R /srv/shiny-server/scaflow/
COPY DB/ /srv/shiny-server/scaflow/DB 

RUN Rscript -e 'install.packages(cran <- c("devtools"))'

RUN Rscript -e 'devtools::install_github("ParkerICI/vite")'RUN Rscript -e 'devtools::install_github(repo = "zunderlab/FLOWMAP")'RUN Rscript -e 'devtools::install_github("DillonHammill/CytoExploreRData")'RUN Rscript -e 'devtools::install_github("DillonHammill/CytoExploreR")'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scaflow

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
