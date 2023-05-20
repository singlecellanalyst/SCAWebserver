FROM rocker/shiny:latest

COPY server.R /srv/shiny-server/scascrna/
  COPY ui.R /srv/shiny-server/scascrna/
  COPY DB/ /srv/shiny-server/scascrna/DB 

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y default-jdk && \
apt-get install -y r-cran-rjava && \
apt-get install -y git libxml2-dev libmagick++-dev gfortran libglpk-dev libgsl-dev libhdf5-dev && \
apt-get install -y tcl tk libglu1 libproj-dev libgdal-dev && \
rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("shiny","BiocManager","devtools","remotes","kml","dplyr","signal","Hmisc","ggplot2","rgl","R.utils"))'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages("https://cran.r-project.org/src/contrib/Archive/clusterCrit/clusterCrit_1.2.8.tar.gz")'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages("https://cran.r-project.org/src/contrib/Archive/akmedoids/akmedoids_1.3.0.tar.gz")'

RUN Rscript -e 'BiocManager::install("multtest",ask = F)'

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(cran <- c("Seurat","dplyr","patchwork","reshape2","reticulate","plotly","viridis","cowplot","plot3D","DDRTree","pheatmap","ggbeeswarm","ggthemes","gplots","gridExtra","ggpubr","lattice","HGNChelper","harmony","ggrepel","grid","mutoss","metap","htmlwidgets","ggnewscale","ggupset","europepmc","dendextend","DT","R.utils","hdf5r","ggridges","scales","xlsx","ggplot2","shiny","shinythemes","shinyWidgets","shinyscreenshot","lme4","terra","ggrastr"))'

RUN Rscript -e 'BiocManager::install(c("scater","ReactomePA","DOSE","enrichplot","multtest","celldex","clusterProfiler","org.Hs.eg.db","scRNAseq","SingleR","ComplexHeatmap","BiocGenerics","DelayedArray","DelayedMatrixStats","limma","S4Vectors","SingleCellExperiment","SummarizedExperiment","batchelor","HDF5Array"),ask = F)'

RUN Rscript -e 'devtools::install_github("cole-trapnell-lab/monocle3",force = T)'
RUN Rscript -e 'devtools::install_github(repo = "satijalab/seurat-wrappers",force = T)'

RUN Rscript -e 'packages <- c("shiny","BiocManager","devtools","remotes","kml","dplyr","signal","Hmisc","ggplot2","rgl","R.utils","Seurat","dplyr","patchwork","reshape2","reticulate","plotly","viridis","cowplot","plot3D","DDRTree","pheatmap","ggbeeswarm","ggthemes","gplots","gridExtra","ggpubr","lattice","HGNChelper","harmony","ggrepel","grid","metap","htmlwidgets","ggnewscale","ggupset","europepmc","dendextend","DT","R.utils","hdf5r","ggridges","scales","xlsx","ggplot2","shiny","shinythemes","shinyWidgets","shinyscreenshot","lme4","terra","ggrastr","scater","ReactomePA","DOSE","enrichplot","multtest","celldex","clusterProfiler","org.Hs.eg.db","scRNAseq","SingleR","ComplexHeatmap","BiocGenerics","DelayedArray","DelayedMatrixStats","limma","S4Vectors","SingleCellExperiment","SummarizedExperiment","batchelor","HDF5Array","monocle3","seurat-wrappers"); lapply(packages, library, character.only = TRUE)'

RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/scascrna

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]