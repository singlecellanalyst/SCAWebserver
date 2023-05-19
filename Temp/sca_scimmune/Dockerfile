FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scascimmune/
COPY ui.R /srv/shiny-server/scascimmune/
COPY DB/ /srv/shiny-server/scascimmune/DB 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y openjdk-17-jre-headless && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("xlsx"))'

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/scascimmune

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
