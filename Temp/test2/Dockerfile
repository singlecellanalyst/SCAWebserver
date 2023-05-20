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
apt-get install -y tcl tk && \
rm -rf /var/lib/apt/lists/*

# apt-get install -y openjdk-17-jre-headless && \

RUN which java
RUN ls -lrt /usr/bin/java
# RUN export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-17-openjdk-amd64/jre
# ENV LD_LIBRARY_PATH=/usr/lib/jvm/java-17-openjdk-amd64/

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(cran <- c("remotes","devtools","BiocManager","R.utils"))'
RUN Rscript -e 'devtools::install_github("thomasp85/patchwork",force = T)'
RUN Rscript -e 'BiocManager::install(c("AnnotationFilter","BiocGenerics","GenomicFeatures","GenomicRanges","IRanges","Rsamtools","S4Vectors","ggbio","motifmatchr","AnnotationDbi","Seurat","GenomeInfoDb","JASPAR2020","TFBSTools","EnsDb.Hsapiens.v86","BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19","ComplexHeatmap"), ask = F)'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(cran <- c("shiny","patchwork","hdf5r","rJava"))'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages("Signac")'

RUN Rscript -e 'packages <- c("remotes","devtools","BiocManager","R.utils","patchwork","AnnotationFilter","BiocGenerics","GenomicFeatures","GenomicRanges","IRanges","Rsamtools","S4Vectors","ggbio","motifmatchr","AnnotationDbi","Seurat","GenomeInfoDb","JASPAR2020","TFBSTools","EnsDb.Hsapiens.v86","BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19","Signac","shiny","patchwork","hdf5r","ComplexHeatmap","rJava"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/test

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]