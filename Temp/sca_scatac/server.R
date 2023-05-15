###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: scATAC Analysis Pipeline
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last-update date: 2023-05-15
# All Rights Reserved
###########################################################################################################################

library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=30000*1024^2)

library("shiny")
library("shinyjs")
library("shinyWidgets")
library("shinyscreenshot")
library("shinythemes")
library("shinyFiles")
library("shinydashboard")
library("shinyalert")
library("webshot")
library("Seurat")
# library("SeuratData")
library("Signac")
library("patchwork")
library("ggpubr")
library("dplyr")
# library("chromVAR")
library("ggplot2")
library("GenomeInfoDb")
library("gridExtra")
library("motifmatchr")
library("JASPAR2020")
library("TFBSTools")
library("cowplot")
library("SeuratWrappers")
library("EnsDb.Hsapiens.v86")
library("BSgenome.Hsapiens.UCSC.hg38")
library("EnsDb.Hsapiens.v75")
library("BSgenome.Hsapiens.UCSC.hg19")
library("hdf5r")
library("ggridges")
library("ggthemes")
library("reshape2")
library("ggrepel")
library("scales")
library("xlsx")
library("qlcMatrix")
library("ggseqlogo")

Sys.setenv("DISPLAY"=":0")
source("DB/SCA_scATAC_RShiny_Functions_V1.0.0.R")
color_conditions <- color_ini()

hpca.se <- readRDS("DB/HumanPrimaryCellAtlasData.RDS")

# example1 <- "DB/SCA_scATAC_Example_From_10X.zip"
example2 <- "DB/SCA_scATAC_Metadata_Example.csv"

ctime <- format(Sys.time(), format = "%Y%m%d%H%M%S", tz = "Europe/Stockholm")

shinyServer(function(input, output, session) {
  
  values <- reactiveValues(proceed = 1)

  observeEvent(input$submit & input$termscheck,{

    if(input$termscheck == FALSE & input$submit){
      showModal(modalDialog("Please agree to our terms and conditions before you proceed", footer=NULL, easyClose = T))
      values$proceed <- 0
    }

    if(input$termscheck == TRUE & input$submit){
      if(values$proceed == 1){
        updateTabsetPanel(session, "nav",selected = "section1")
      }else{
        values$proceed <- 1
      }
    }
  })

  # output$downloadExample1 <- downloadHandler(
  #   filename = "SCA_scATAC_Example_From_10X.zip",
  #   content = function(file) {
  #     file.copy(example1, file)
  #   },
  #   contentType = "application/zip"
  # )

  output$downloadExample2 <- downloadHandler(
    filename = "SCA_scATAC_Metadata_Example.csv",
    content = function(file) {
      file.copy(example2, file)
    },
    contentType = "text/csv"
  )

  data <- reactive({
    req(input$files)
    showModal(modalDialog("Processing..depending on the number of samples as well as file sizes, time of processing may varies. Please wait patiently for results to be delivered.", footer=NULL, easyClose = F))
    # results <- readRDS("DB/RESULTS_EXAMPLE.RDS")
    print(input$phenodata$datapath)
    pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "scATAC", isDir = T)
    pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, sep = "_")
    sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SID)))
    names(sample_colors) <- unique(pheno_data$SID)

    print(input$files)
    print(input$files$datapath)
    project_name <- ifelse(input$project_name == "","My_Project",gsub("\\s+|\\(|\\)|-|\\/|\\?","",input$project_name))
    print(project_name)
    cdir <- paste(getwd(),"/",project_name,ctime, "/",sep = "")
    dir.create(cdir)

    sample_files <- unzip(input$files$datapath, exdir = cdir)
    sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]

    removeModal()

    #######################################################################################################################################
    showModal(modalDialog("Preparing..", footer=NULL, easyClose = F))

    pheno_peak <- list.files(cdir, pattern = ".*peak.*h5", full.names = T, ignore.case = T, recursive = T)
    pheno_singlecell <- list.files(cdir, pattern = ".*singlecell\\.csv$|.*per_barcode_metrics\\.csv$", full.names = T, ignore.case = T, recursive = T)
    pheno_fragments <- list.files(cdir, pattern = ".*fragment.*tsv.*", full.names = T, ignore.case = T, recursive = T)
    pheno_fragments <- pheno_fragments[grep("\\.tbi", pheno_fragments, ignore.case = T, invert = T)]
    error_files <- list.files(cdir, pattern = "gz\\.tsv$", full.names = T, ignore.case = T)

    if(length(error_files) > 0){
      for(i in 1:length(error_files)){
        current_corrected_file <- gsub("gz\\.tsv","gz",error_files[i])
        system(paste("mv ", error_files[i], " ",current_corrected_file, sep = ""))
      }
    }

    data <- NULL
    sample_names <- NULL
    results <- NULL

    removeModal()

    for(i in 1:nrow(pheno_data)){
      showModal(modalDialog(paste("Running ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      sample_names <- c(sample_names, pheno_data[i,"SID"])
      current <- NULL
      current <- pheno_peak[grep(pheno_data[i,"FILE"], pheno_peak, ignore.case = T)]
      current <- Read10X_h5(current)
      current_singlecell <- NULL
      current_singlecell <- pheno_singlecell[grep(pheno_data[i,"FILE"], pheno_singlecell, ignore.case = T)]

      if(length(current_singlecell) > 0){
        current_meta <- "YES"
        current_singlecell <- read.csv(file = current_singlecell, header = TRUE, row.names = 1)
      }else{
        current_meta <- "NO"
      }

      if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
        ref_genome <- "hg19"
        load("DB/hg19_EnsDb.Hsapiens.v75.RData")
        # load("/mnt/REFERENCE/hg19_EnsDb.Hsapiens.v75.RData")
        seqlevelsStyle(ref_annot) <- "UCSC"
        genome(ref_annot) <- "hg19"
      }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
        ref_genome <- "GRCh38.p12"
        load("DB/hg38_EnsDb.Hsapiens.v86.RData")
        # load("/mnt/REFERENCE/hg38_EnsDb.Hsapiens.v86.RData")
        seqlevelsStyle(ref_annot) <- "UCSC"
        genome(ref_annot) <- "GRCh38.p12"
      }else{
        stop("You are supplying non-human samples. Please try again.")
      }

      removeModal()
      showModal(modalDialog(paste("Constructing chromatin assay of ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))

      current_assay <- CreateChromatinAssay(
        counts = current,
        sep = c(unique(gsub(".*?([-.:])[0-9]+([-.:])[0-9]+","\\1",row.names(current))),
                unique(gsub(".*?([-.:])[0-9]+([-.:])[0-9]+","\\2",row.names(current)))),
        genome = ref_genome,
        fragments = pheno_fragments[grep(pheno_data[i,"FILE"], pheno_fragments, ignore.case = T)],
        min.cells = 10,
        min.features = 200)

      if(current_meta == "YES"){
        current_seurat <- CreateSeuratObject(counts = current_assay, assay = "peaks",meta.data = current_singlecell)
        rm(current_singlecell)
      }else{
        current_seurat <- CreateSeuratObject(counts = current_assay, assay = "peaks")
      }

      current_seurat <- add_names(current_seurat, pheno_data[i,"SID"], current_ident = "orig.ident")

      rm(current_assay)
      rm(current)

      removeModal()
      showModal(modalDialog(paste("Finding nucelosome signal for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))

      DefaultAssay(current_seurat) <- 'peaks'
      Annotation(current_seurat) <- ref_annot
      current_seurat <- NucleosomeSignal(object = current_seurat)
      removeModal()
      showModal(modalDialog(paste("Calculating TSS enrichment for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      current_seurat <- TSSEnrichment(object = current_seurat, fast = FALSE)
      removeModal()
      showModal(modalDialog(paste("Calculating fragments for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      if(length(grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat$pct_reads_in_peaks <- current_seurat@meta.data[,grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)] / current_seurat@meta.data[,grep("passed_filters", colnames(current_seurat@meta.data),ignore.case = T)] * 100
        if(length(grep("blacklist_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
          current_seurat$blacklist_ratio <- current_seurat@meta.data[,grep("blacklist_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)] / current_seurat@meta.data[,grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)]
        }
      }
      current_seurat$high.tss <- ifelse(current_seurat@meta.data[,grep("TSS.enrichment",colnames(current_seurat@meta.data), ignore.case = T)] > 2, 'High', 'Low')

      p <- NULL
      p <- TSSPlot(current_seurat, group.by = 'high.tss') + NoLegend()+ scale_color_manual(values = color_conditions$alternate, length(current_seurat$high.tss))

      results[['p1plots']][[i]] <- adjust_theme(p, legend = "none")
      names(results[['p1plots']])[i] <- pheno_data[i,"SID"]

      p <- NULL
      current_col <- grep("nucleosome_signal",colnames(current_seurat@meta.data), ignore.case = T)
      current_seurat$nucleosome_group <- ifelse(current_seurat@meta.data[,current_col] > 4, 'NS > 4', 'NS < 4')
      current_groups <- unique(current_seurat$nucleosome_group)
      if(length(current_groups) == 1){
        p <- FragmentHistogram(object = current_seurat) +
          scale_fill_manual(values = color_conditions$manycolors, length(current_seurat$nucleosome_group))+ggtitle(paste("FRAGMENT LENGTH PERIODICITY: ", current_groups, sep = ""))
      }else{
        p <- FragmentHistogram(object = current_seurat, group.by = 'nucleosome_group') +
          scale_fill_manual(values = color_conditions$manycolors, length(current_seurat$nucleosome_group))+ggtitle("FRAGMENT LENGTH PERIODICITY")
      }

      results[['p2plots']][[i]] <- adjust_theme(p, legend = "none")
      names(results[['p2plots']])[i] <- pheno_data[i,"SID"]

      colnames(current_seurat@meta.data)[grep("pct_reads_in_peaks", colnames(current_seurat@meta.data), ignore.case = T)] <- "PCT_Reads_in_Peaks"
      colnames(current_seurat@meta.data)[grep("blacklist_ratio", colnames(current_seurat@meta.data), ignore.case = T)] <- "Blacklist_Ratio"
      colnames(current_seurat@meta.data)[grep("peak_region_fragments", colnames(current_seurat@meta.data), ignore.case = T)] <- "Peak_Region_Fragments"
      colnames(current_seurat@meta.data)[grep("TSS.enrichment", colnames(current_seurat@meta.data), ignore.case = T)] <- "TSS_Enrichment"
      colnames(current_seurat@meta.data)[grep("nucleosome_signal", colnames(current_seurat@meta.data), ignore.case = T)] <- "Nucleosome_Signal"

      current_features <- colnames(current_seurat@meta.data)[grep("PCT_Reads_in_Peaks|Peak_Region_Fragments$|TSS_Enrichment$|Blacklist_Ratio$|Nucleosome_Signal$", colnames(current_seurat@meta.data), ignore.case = T)]

      p <- NULL

      results[['p3plots']][[i]] <- violin_plot(current_seurat, features = current_features, ncol = length(current_features),
                                               col = sample_colors, x_lab = "SID", log_status = TRUE)
      names(results[['p3plots']])[i] <- pheno_data[i,"SID"]

      if(length(grep("Peak_Region_Fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = Peak_Region_Fragments > 3000 & Peak_Region_Fragments < 20000)
      }

      if(length(grep("^nCount_peaks$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = nCount_peaks < 100000 & nCount_peaks > 1000)
      }

      if(length(grep("^PCT_Reads_in_Peaks$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = PCT_Reads_in_Peaks > 15)
      }

      if(length(grep("^Blacklist_Ratio$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = Blacklist_Ratio < 0.05)
      }

      if(length(grep("^Nucleosome_Signal$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = Nucleosome_Signal < 4)
      }

      if(length(grep("^TSS_Enrichment$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = TSS_Enrichment > 2)
      }

      if(length(grep("^nCount_RNA$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = nCount_RNA < 25000 & nCount_RNA > 1000)
      }

      if(length(grep("^nFeature_RNA$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
      }

      if(length(grep("^Percent_Mito$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat <- subset(x = current_seurat, subset = Percent_Mito < 5)
      }

      removeModal()
      showModal(modalDialog(paste("Running dimension reduction and clustering for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      DefaultAssay(current_seurat) <- "peaks"
      current_seurat <- RunTFIDF(current_seurat)
      current_seurat <- FindTopFeatures(current_seurat, min.cutoff = 'q0')
      current_seurat <- RunSVD(current_seurat)
      results[['p4plots']][[i]] <- DepthCor(current_seurat) + theme_bw(base_size = 20)
      names(results[['p4plots']])[i] <- pheno_data[i,"SID"]

      dim_i <- 2 # Skip first component with high correlation
      m <- ifelse((length(current_seurat@assays$peaks@var.features)-1) < 30, (length(current_seurat@assays$peaks@var.features) - 1), 30)
      current_seurat <- RunUMAP(object = current_seurat, reduction = 'lsi', dims = dim_i:m)
      current_seurat <- FindNeighbors(object = current_seurat, reduction = 'lsi', dims = dim_i:m)
      current_seurat <- FindClusters(object = current_seurat, verbose = FALSE, algorithm = 3)

      plotx <- NULL
      plotx <- gen10x_plotx(current_seurat, selected = c("UMAP"), include_meta = T)
      plotx$CLUSTER <- plotx$seurat_clusters
      cluster_colors <- gen_colors(color_conditions$tenx, length(unique(plotx$CLUSTER)))
      names(cluster_colors) <- sort(unique(plotx$CLUSTER), decreasing = F)

      p <- NULL
      p <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "CLUSTER", plot_title = project_name, col = cluster_colors, annot = TRUE,
                        numeric = T, legend_position = "right", point_size = 1, label_size = 8)
      results[['p5plots']][[i]] <- p
      names(results[['p5plots']])[i] <- pheno_data[i,"SID"]

      removeModal()
      showModal(modalDialog(paste("Finding top peaks for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      p_thresh <- 0.05
      DefaultAssay(current_seurat) <- "peaks"
      da_peaks <- NULL
      da_peaks <- FindAllMarkers(
        object = current_seurat,
        min.pct = 0.5,
        test.use = 'LR',
        return.thresh = p_thresh,
        latent.vars = if(is.null(current_seurat@meta.data$Peak_Region_Fragments)) newobj <- NULL else newobj <- "Peak_Region_Fragments")

      if(nrow(da_peaks) > 0){
        da_peaks <- da_peaks[da_peaks$p_val_adj < p_thresh,]
        wt <- colnames(da_peaks)[grep("log.*FC", colnames(da_peaks), ignore.case = T)]

        da_peaks <- da_peaks[order(da_peaks[,wt], decreasing = TRUE),]
        results[['p6data']][[i]] <- da_peaks[,c(grep("gene", colnames(da_peaks), ignore.case = T),
                                                grep("gene", colnames(da_peaks), ignore.case = T, invert = T))]
        names(results[['p6data']])[i] <- pheno_data[i,"SID"]

        n <- 6
        topn <- split(da_peaks, da_peaks$cluster)
        topn <- lapply(topn, function(x){
          x <- x[order(x[,wt]),]
          x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
        })
        topn <- do.call(rbind.data.frame, topn)

        closest_genes <- NULL
        cclusters <- unique(topn$cluster)
        for(k in 1:length(cclusters)){
          current <- unique(topn[which(topn$cluster == cclusters[k]),"gene"])
          current <- data.frame(cluster = cclusters[k], ClosestFeature(current_seurat, regions = current))
          closest_genes <- rbind(closest_genes, current)
        }
        results[['p8data']][[i]] <- closest_genes
        names(results[['p8data']])[i] <- pheno_data[i,"SID"]
      }
      removeModal()
      showModal(modalDialog(paste("Done with top peaks search for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      removeModal()
      showModal(modalDialog(paste("Deciphering gene activities for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      gene_activities <- GeneActivity(current_seurat)
      current_seurat[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
      DefaultAssay(current_seurat) <- "ACTIVITY"
      current_seurat <- NormalizeData(current_seurat)
      current_seurat <- ScaleData(current_seurat, features = rownames(current_seurat))
      Idents(current_seurat) <- "seurat_clusters"
      removeModal()
      showModal(modalDialog(paste("Finding top genes for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      current_data_markers <- FindAllMarkers(current_seurat, min.pct = 0.5, logfc.threshold = 0.25)
      p_thresh <- 0.05
      current_data_markers <- current_data_markers[current_data_markers$p_val_adj < p_thresh,]
      current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]
      results[['p9data']][[i]] <- current_data_markers
      names(results[['p9data']])[i] <- pheno_data[i,"SID"]

      DefaultAssay(current_seurat) <- "ACTIVITY"
      wt <- colnames(da_peaks)[grep("log.*FC", colnames(da_peaks), ignore.case = T)]
      top1 <- split(da_peaks, da_peaks$cluster)
      top1 <- lapply(top1, function(x){
        x <- x[order(x[,wt]),]
        if(nrow(x) > 0){
          x <- x[1:ifelse(nrow(x)>1, 1, nrow(x)),]
        }
      })
      top1 <- do.call(rbind.data.frame, top1)

      top1_derivedrna <- split(current_data_markers, current_data_markers$cluster)
      top1_derivedrna <- lapply(top1_derivedrna, function(x){
        x <- x[order(x[,wt]),]
        if(nrow(x) > 0){
          x <- x[1:ifelse(nrow(x)>1, 1, nrow(x)),]
        }
      })

      top1_derivedrna <- do.call(rbind.data.frame, top1_derivedrna)

      removeModal()
      showModal(modalDialog(paste("Generating regional statistics for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      DefaultAssay(current_seurat) <- "peaks"
      main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
      keep.peaks <- which(as.character(seqnames(granges(current_seurat))) %in% main.chroms)
      current_seurat[["peaks"]] <- subset(current_seurat[["peaks"]], features = rownames(current_seurat[["peaks"]])[keep.peaks])
      rm(keep.peaks)
      rm(main.chroms)
      if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
        current_seurat <- RegionStats(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg19)
      }else if(length(grep("hg38|grch38|38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
        current_seurat <- RegionStats(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg38)
      }

      removeModal()
      showModal(modalDialog(paste("Linking peaks to genes for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      current_seurat <- LinkPeaks(
        object = current_seurat,
        peak.assay = "peaks",
        expression.assay = "ACTIVITY",
        genes.use = unlist(top1_derivedrna$gene))

      clusters <- sort(as.numeric(as.character(unique(top1$cluster))),decreasing = F)
      DefaultAssay(current_seurat) <- "peaks"
      for(k in 1:length(clusters)){
        if(length(which(unlist(top1$cluster) == clusters[k])) > 0){
          results[['p10plots']][[length(results[['p10plots']])+1]] <- CoveragePlot(
            object = current_seurat,
            region = unlist(top1[which(top1$cluster == clusters[k]),"gene"]),
            features = unlist(top1_derivedrna[which(top1_derivedrna$cluster == clusters[k]),"gene"]),
            expression.assay = "ACTIVITY",
            extend.upstream = 10000,
            extend.downstream = 10000)+
            plot_annotation(title = paste("TOP PEAK AND GENE ACTIVITY\n",pheno_data[i,"SID"], ": CLUSTER ", clusters[k], sep = ""),
                            theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
          names(results[['p10plots']])[length(results[['p10plots']])] <- paste(pheno_data[i,"SID"],"|","CLUSTER_",clusters[k],sep = "")
        }
      }

      n <- 6
      topn_genes <- split(current_data_markers, current_data_markers$cluster)
      topn_genes <- lapply(topn_genes, function(x){
        x <- x[order(x[,wt]),]
        x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
      })
      topn_genes <- do.call(rbind.data.frame, topn_genes)
      # results[['p11data']][[i]] <- topn_genes
      # names(results[['p11data']])[i] <- pheno_data[i,"SID"]

      removeModal()
      showModal(modalDialog(paste("Add motif information for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      current_assay <- "peaks"
      used_type <- "ACTIVITY"
      DefaultAssay(current_seurat) <- "peaks"
      chosen_genes <- NULL
      position_freq <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
      if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
        current_seurat <- AddMotifs(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = position_freq)
        motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
        if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
          chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
          chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
          current_seurat <- Footprint(object = current_seurat, motif.name = chosen_genes,
                                      genome = BSgenome.Hsapiens.UCSC.hg19)
        }
      }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
        current_seurat <- AddMotifs(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = position_freq)
        motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
        if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
          chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
          chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
          current_seurat <- Footprint(object = current_seurat, motif.name = chosen_genes,
                                      genome = BSgenome.Hsapiens.UCSC.hg38)
        }
      }

      p <- NULL
      p <- PlotFootprint(current_seurat, features = chosen_genes, group.by	= "seurat_clusters")
      p <- p + plot_layout(ncol = 2)+
        plot_annotation(title = paste(pheno_data[i,"SID"], ": TOP ",used_type," GENES \n(Top group with highest accessibility in motif flanking region are labeled)", sep = ""),
                        theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
      results[['p12plots']][[i]] <- p
      names(results[['p12plots']])[i] <- pheno_data[i,"SID"]

      p_thresh <- 0.05
      DefaultAssay(current_seurat) <- "peaks"

      enriched.motifs <- FindMotifs(object = current_seurat,features = da_peaks[which(da_peaks$p_val_adj < p_thresh),"gene"])
      results[['p13data']][[i]] <- enriched.motifs
      names(results[['p13data']])[i] <- pheno_data[i,"SID"]

      n <- 10
      results[['p14plots']][[i]] <- MotifPlot(object = current_seurat, motifs = rownames(enriched.motifs)[1:n])+
        plot_annotation(title = paste(pheno_data[i,"SID"], "\nPOSITION WEIGHT MATRICES OF TOP ENRICHED MOTIFS", sep = ""),
                        theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
      results[['p14plots']][[i]] <- adjust_theme(results[['p14plots']][[i]], xsize = 12, title_size = 20, strip_size = 15)
      names(results[['p14plots']])[i] <- pheno_data[i,"SID"]

      removeModal()
      showModal(modalDialog(paste("Computing peak set variability for ", pheno_data[i,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
      # Per-cell motif activity score
      # if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      #   current_seurat <- RunChromVAR(object = current_seurat,genome = BSgenome.Hsapiens.UCSC.hg19)
      # }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      #   current_seurat <- RunChromVAR(object = current_seurat,genome = BSgenome.Hsapiens.UCSC.hg38)
      # }
      # 
      # n <- 6
      # DefaultAssay(current_seurat) <- 'chromvar'
      # p <- NULL
      # p <- FeaturePlot(
      #   object = current_seurat,
      #   features = enriched.motifs$motif[1:n],
      #   min.cutoff = 'q10',
      #   max.cutoff = 'q90',
      #   pt.size = 0.1, label = T,
      #   cols = c("green", "blue")) +
      #   plot_annotation(title = "Motif Activities",
      #                   theme = theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold")))
      # 
      # results[['p15plots']][[i]] <- (results[['p5plots']][[i]] / p) +
      #   plot_annotation(title = paste(pheno_data[i,"SID"], "\nPeaks VS Top Enriched Motifs", sep = ""),
      #                   theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
      # names(results[['p15plots']])[i] <- pheno_data[i,"SID"]

      DefaultAssay(current_seurat) <- "peaks"
      n <- 1000
      topn <- split(da_peaks, da_peaks$cluster)
      topn <- lapply(topn, function(x){
        x <- x[order(x[,wt]),]
        x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
      })
      topn <- do.call(rbind.data.frame, topn)
      n <- 5000
      topn <- unique(unlist(topn$gene)[1:ifelse(nrow(topn) < n, nrow(topn), n)])
      # current_seurat <- subset(current_seurat, features = topn)
      results[['data']][[i]] <- current_seurat
      names(results[['data']])[i] <- pheno_data[i,"SID"]
      rm(current_seurat)
      rm(closest_genes)
      rm(position_freq)
      # rm(top_FC)
      # rm(top.da.peak)
      removeModal()
    }

    data_ref <- unique(toupper(pheno_data[,"REF_GENOME"]))
    if((length(results[['data']]) > 1) & (length(data_ref) == 1)){
      showModal(modalDialog(paste("Merging samples..", sep = ""), footer=NULL, easyClose = F))
      combined_peaks <- NULL
      for(i in 1:(length(results[['data']]) - 1)){
        combined_peaks <- reduce(x = c(results[['data']][[i]]@assays$peaks@ranges, results[['data']][[i+1]]@assays$peaks@ranges))
      }

      peakwidths <- width(combined_peaks)
      combined_peaks <- combined_peaks[peakwidths  < 10000 & peakwidths > 20]

      for(i in 1:length(results[['data']])){

        current <- FeatureMatrix(
          fragments = results[['data']][[i]]@assays$peaks@fragments,
          features = combined_peaks,
          cells = colnames(results[['data']][[i]]))

        current <- CreateChromatinAssay(current, fragments = results[['data']][[i]]@assays$peaks@fragments)
        results[['data']][[i]] <- CreateSeuratObject(current, assay = "peaks")
        results[['data']][[i]]$Sample <- sample_names[i]

      }

      rm(current)

      data <- merge(x = results[['data']][[1]], y = results[['data']][2:length(results[['data']])],
                    add.cell.ids = sample_names)
      removeModal()
      showModal(modalDialog(paste("Running dimension reduction and clustering for merged samples..", sep = ""), footer=NULL, easyClose = F))

      n <- ifelse((length(data@assays$peaks@var.features)-1) < 50, (length(data@assays$peaks@var.features) - 1), 50)
      data <- RunTFIDF(data)
      data <- FindTopFeatures(data, min.cutoff = 20)
      data <- RunSVD(data) # , n = n
      data <- RunUMAP(data, reduction = 'lsi', dims = 2:n)
      # data <- RunTSNE(data, reduction = 'lsi', check_duplicates = FALSE) # dims = 2:n,
      data <- FindNeighbors(data, reduction = 'lsi') # dims = 2:n,
      data <- FindClusters(data, resolution = 0.8, n.iter = 1000)

      removeModal()
      showModal(modalDialog(paste("Finding top peaks for merged samples..", sep = ""), footer=NULL, easyClose = F))
      data_markers <- FindAllMarkers(data, min.pct = 0.5, logfc.threshold = 0.25)

      p_thresh <- 0.05
      data_markers <- data_markers[data_markers$p_val_adj < p_thresh,]
      data_markers <- data_markers[order(data_markers$p_val_adj, decreasing = F),]
      results[['data_markers']] <- data_markers

      plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                          UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                          Sample = data$Sample, CLUSTER = data$seurat_clusters)
      plotx$Sample <- factor(plotx$Sample)
      plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(unique(as.numeric(as.character(plotx$CLUSTER)))))

      cluster_colors <- gen_colors(color_conditions$tenx, length(unique(plotx$CLUSTER)))
      names(cluster_colors) <- sort(unique(plotx$CLUSTER), decreasing = F)

      p1_umap <- NULL
      p1_umap <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "Sample", plot_title = "All Samples Integrated scATAC UMAP\nBy Samples", col = sample_colors, annot = TRUE, legend_position = "right", point_size = 1, numeric = F)

      p2_umap <- NULL
      p2_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "CLUSTER", plot_title = "All Samples Integrated scATAC UMAP\nBy Clusters", col = cluster_colors, annot = TRUE, legend_position = "right", point_size = 1, numeric = T)

      results[['p16plots']] <- p1_umap+p2_umap

      wt <- colnames(results[['data_markers']])[grep("log.*FC", colnames(results[['data_markers']]), ignore.case = T)]
      top1 <- results[['data_markers']] %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
      top1 <- top1[order(top1$cluster, decreasing = F),]

      plots <- NULL
      for(k in 1:length(top1$gene)){

        plots[[k]] <- FeaturePlot(data, features = top1$gene[k], cols = c("green", "blue"),
                                  label = T, pt.size = 0.1, max.cutoff = 'q95')+
          ggtitle(paste("CLUSTER: ", top1$cluster[k], "\n",top1$gene[k],sep = ""))
      }

      results[['p17plots']] <- plots

      n <- 10
      wt <- colnames(results[['data_markers']])[grep("log.*FC", colnames(results[['data_markers']]), ignore.case = T)]
      top_n <- results[['data_markers']] %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
      current_clusters <- as.numeric(as.character(sort(unique(results[['data_markers']]$cluster))))
      for(k in 1:length(current_clusters)){
        if(length(which(unlist(top_n$cluster) == current_clusters[k])) > 0){
          current <- top_n[which(top_n$cluster == current_clusters[k]),]
          current <- current[order(current$p_val_adj, decreasing = F),]
          results[['p18plots']][[length(results[['p18plots']])+1]] <- VlnPlot(data, features = current$gene, pt.size = 0,
                                                                              # log = TRUE,
                                                                              cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) +
            plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", current_clusters[k], sep = ""),
                            theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
          names(results[['p18plots']])[length(results[['p18plots']])] <- paste("CLUSTER_",current_clusters[k], sep = "")
        }
      }

      if(length(grep("hg19",data_ref, ignore.case = T)) > 0){
        ref_genome <- "hg19"
        load("DB/hg19_EnsDb.Hsapiens.v75.RData")
        seqlevelsStyle(ref_annot) <- "UCSC"
        genome(ref_annot) <- "hg19"
      }else if(length(grep("hg38|grch38",data_ref, ignore.case = T)) > 0){
        ref_genome <- "GRCh38.p12"
        load("DB/hg38_EnsDb.Hsapiens.v86.RData")
        seqlevelsStyle(ref_annot) <- "UCSC"
        genome(ref_annot) <- "GRCh38.p12"
      }

      Annotation(data) <- ref_annot
      rm(ref_annot)
      removeModal()
      showModal(modalDialog(paste("Calculating gene activities for merged samples..", sep = ""), footer=NULL, easyClose = F))
      gene_activities <- GeneActivity(data)
      data[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
      rm(gene_activities)

      DefaultAssay(data) <- "ACTIVITY"

      removeModal()
      showModal(modalDialog(paste("Finding top genes for merged samples..", sep = ""), footer=NULL, easyClose = F))
      data_activity_markers <- FindAllMarkers(data, min.pct = 0.5, logfc.threshold = 0.25)
      p_thresh <- 0.05
      data_activity_markers <- data_activity_markers[data_activity_markers$p_val_adj < p_thresh,]
      data_activity_markers <- data_activity_markers[order(data_activity_markers$p_val_adj, decreasing = F),]
      results[['data_activity_markers']] <- data_activity_markers

      wt <- colnames(data_activity_markers)[grep("log.*FC", colnames(data_activity_markers), ignore.case = T)]
      top1_activities <- data_activity_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
      top1_activities <- top1_activities[order(top1_activities$cluster, decreasing = F),]

      top1 <- data_activity_markers %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
      top1 <- top1[order(top1$cluster, decreasing = F),]

      clusters <- unique(sort(as.numeric(as.character(top1$cluster))))

      DefaultAssay(data) <- "peaks"
      for(k in 1:length(clusters)){
        if(length(which(unlist(top1$cluster) == clusters[k])) > 0){
          results[['p19plots']][[length(results[['p19plots']])+1]] <- CoveragePlot(
            object = data,
            expression.assay = "ACTIVITY",
            region = unlist(top1[which(top1$cluster == clusters[k]),"gene"]),
            features = unlist(top1_activities[which(top1_activities$cluster == clusters[k]),"gene"]),
            extend.upstream = 10000,
            extend.downstream = 10000)+
            # scale_color_manual(values = gen_colors(color_conditions$tenx, length(unique(data$Sample))))+
            plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", clusters[k], sep = ""),
                            theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
          names(results[['p19plots']])[length(results[['p19plots']])] <- paste("CLUSTER_",clusters[k], sep = "")
        }
      }

      removeModal()
    }

    #######################################################################################################################################
    # saveRDS(results,"results.RDS")
    showModal(modalDialog("Preparing plots..", footer=NULL, easyClose =  F))
    updateSelectInput(session, inputId = 'p1id', label = 'Choose a sample to display', choices = names(results[['p1plots']]), selected = names(results[['p1plots']])[1])
    updateSelectInput(session, inputId = 'p2id', label = 'Choose a sample to display', choices = names(results[['p2plots']]), selected = names(results[['p2plots']])[1])
    updateSelectInput(session, inputId = 'p3id', label = 'Choose a sample to display', choices = names(results[['p3plots']]), selected = names(results[['p3plots']])[1])
    updateSelectInput(session, inputId = 'p4id', label = 'Choose a sample to display', choices = names(results[['p4plots']]), selected = names(results[['p4plots']])[1])
    updateSelectInput(session, inputId = 'p5id', label = 'Choose a sample to display', choices = names(results[['p5plots']]), selected = names(results[['p5plots']])[1])
    updateSelectInput(session, inputId = 'p6id', label = 'Choose a sample to display', choices = names(results[['p6data']]), selected = names(results[['p6data']])[1])
    # updateSelectInput(session, inputId = 'p7id', label = 'Choose a sample to display', choices = names(results[['p7data']]), selected = names(results[['p7data']])[1])
    updateSelectInput(session, inputId = 'p8id', label = 'Choose a sample to display', choices = names(results[['p8data']]), selected = names(results[['p8data']])[1])
    updateSelectInput(session, inputId = 'p9id', label = 'Choose a sample to display', choices = names(results[['p9data']]), selected = names(results[['p9data']])[1])
    p10samples <- gsub("(.*)\\|.*","\\1",names(results[['p10plots']]))
    p10clusters <- sort(as.numeric(as.character(unique(gsub(".*\\|CLUSTER_(.*)","\\1",names(results[['p10plots']]))))), decreasing = F)
    updateSelectInput(session, inputId = 'p10id', label = 'Choose a sample to display', choices = p10samples, selected = p10samples[1])
    updateSliderInput(session, inputId ="p11id", label = 'Choose a cluster to display', value = 0, min = min(p10clusters), max = max(p10clusters), step = 1)
    updateSelectInput(session, inputId = 'p12id', label = 'Choose a sample to display', choices = names(results[['p12plots']]), selected = names(results[['p12plots']])[1])
    updateSelectInput(session, inputId = 'p13id', label = 'Choose a sample to display', choices = names(results[['p13data']]), selected = names(results[['p13data']])[1])
    updateSelectInput(session, inputId = 'p14id', label = 'Choose a sample to display', choices = names(results[['p14plots']]), selected = names(results[['p14plots']])[1])
    # updateSelectInput(session, inputId = 'p15id', label = 'Choose a sample to display', choices = names(results[['p15plots']]), selected = names(results[['p15plots']])[1])
    p18clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)","\\1",names(results[['p18plots']]))))), decreasing = F)
    updateSliderInput(session, inputId ="p18id", label = 'Choose a cluster to display', value = 0, min = min(p18clusters), max = max(p18clusters), step = 1)
    p19clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)","\\1",names(results[['p19plots']]))))), decreasing = F)
    updateSliderInput(session, inputId ="p19id", label = 'Choose a cluster to display', value = 0, min = min(p19clusters), max = max(p19clusters), step = 1)

    print("Completed!")
    removeModal()
    return(results)
  })

  output$p1plot <- renderPlot({ #renderPlotly
    showModal(modalDialog("Plotting figure 1..", footer=NULL, easyClose = T))
    p <- data()[['p1plots']][[input$p1id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)

  observeEvent(input$p1plot, {
    screenshot(id="p1plot", filename = paste("1SCA_scATAC_TSS_ENRICHMENT_SCORES_",input$p1id,"_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p2plot <- renderPlot({
    showModal(modalDialog("Plotting figure 2..", footer=NULL, easyClose = T))
    p <- data()[['p2plots']][[input$p2id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)

  observeEvent(input$p2plot, {
    screenshot(id="p2plot", filename = paste("2SCA_scATAC_FRAGMENT_LENGTH_PERIODICITY_",input$p2id,"_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p3plot <- renderPlot({
    showModal(modalDialog("Plotting figure 3..", footer=NULL, easyClose =  T))
    p <- data()[['p3plots']][[input$p3id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)

  observeEvent(input$p3plot, {
    screenshot(id="p3plot", filename = paste("3SCA_scATAC_QUALITY_CONTROL_PEAKS_LOGSCALE_",input$p3id,"_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p4plot <- renderPlot({
    showModal(modalDialog("Plotting figure 4..", footer=NULL, easyClose =  T))
    p <- data()[['p4plots']][[input$p4id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)

  observeEvent(input$p4plot, {
    screenshot(id="p4plot", filename = paste("4SCA_scATAC_CORRELATION_SEQUENCING_DEPTH_LSI_COMPONENTS_",input$p4id,"_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p5plot <- renderPlot({
    showModal(modalDialog("Plotting figure 5..", footer=NULL, easyClose = T))
    p <- data()[['p5plots']][[input$p5id]]
    removeModal()
    return(p)
  }, height = 600, width = 900)

  observeEvent(input$p5plot, {
    screenshot(id="p5plot", filename = paste("5SCA_scATAC_UMAP_scATAC_AUTOCLUST_",input$p5id,"_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p6table <- renderDT(data()[['p6data']][[input$p6id]],
                             filter = "top",
                             style="bootstrap",
                             rownames = F,
                             options = list(pageLength = 10))

  output$p8table <- renderDT(data()[['p8data']][[input$p8id]],
                             filter = "top",
                             style="bootstrap",
                             rownames = F,
                             options = list(pageLength = 10))

  output$p9table <- renderDT(data()[['p9data']][[input$p9id]],
                             filter = "top",
                             style="bootstrap",
                             rownames = F,
                             options = list(pageLength = 10))

  output$p10plot <- renderPlot({
    showModal(modalDialog("Plotting figure 10..", footer=NULL, easyClose = T))
    p <- data()[['p10plots']][[paste(input$p10id,"|CLUSTER_",input$p11id, sep = "")]]
    removeModal()
    return(p)
  }, height = 700, width = 900)

  observeEvent(input$p10plot, {
    screenshot(id="p10plot", filename = paste("10SCA_scATAC_ATAC_Tn5_INSERTION_COVERAGE_PEAK_TOP_ACTIVITY_GENE_",input$p10id,"_CLUSTER_",input$p11id, sep = ""), scale = 2)
  })

  output$p12plot <- renderPlot({
    showModal(modalDialog("Plotting figure 12..", footer=NULL, easyClose = T))
    p <- data()[['p12plots']][[input$p12id]]
    removeModal()
    return(p)
  }, height = 900, width = 900)

  observeEvent(input$p12plot, {
    screenshot(id="p12plot", filename = paste("12SCA_scATAC_PEAKS_TOP_ACTIVITY_GENES_FOOTPRINTING_",input$p12id,"_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p13table <- renderDT(data()[['p13data']][[input$p13id]],
                              filter = "top",
                              style="bootstrap",
                              rownames = F,
                              options = list(pageLength = 10))

  output$p14plot <- renderPlot({
    showModal(modalDialog("Plotting figure 14..", footer=NULL, easyClose = T))
    p <- data()[['p14plots']][[input$p14id]]
    removeModal()
    return(p)
  }, height = 700, width = 900)

  observeEvent(input$p14plot, {
    screenshot(id="p14plot", filename = paste("14SCA_scATAC_TOP_GENES_ENRICHED_MOTIFS_",input$p14id, sep = ""), scale = 2)
  })

  # output$p15plot <- renderPlot({
  #   showModal(modalDialog("Plotting figure 15..", footer=NULL, easyClose = T))
  #   p <- data()[['p15plots']][[input$p15id]]
  #   removeModal()
  #   return(p)
  # }, height = 1400, width = 900)

  # observeEvent(input$p15plot, {
  #   screenshot(id="p15plot", filename = paste("15SCA_scATAC_UMAP_VS_TOP_ENRICHED_MOTIFS_",input$p15id, sep = ""), scale = 2)
  # })

  output$p16plot <- renderPlot({
    if(length(grep("ggplot", class(data()[['p16plots']]), ignore.case = T)) > 0){
      showModal(modalDialog("Plotting figure 16..", footer=NULL, easyClose = T))
      p <- data()[['p16plots']]
      removeModal()
    }else{
      text <- paste("\n   No multiple samples\n",
                    "       are provided for merging in this project.")
      p <- ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }
    return(p)
  }, height = 450, width = 900)

  observeEvent(input$p16plot, {
    screenshot(id="p16plot", filename = paste("16SCA_scATAC_UMAP_INTEGRATED_ATAC_SAMPLES_AUTOCLUST_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p17plot <- renderPlot({
    if(length(grep("ggplot", class(data()[['p17plots']]), ignore.case = T)) > 0){
      showModal(modalDialog("Plotting figure 17..", footer=NULL, easyClose = T))
      p <- data()[['p17plots']]
      removeModal()
      print(do.call("grid.arrange", c(p, ncol=4)))
    }else{
      text <- paste("\n   No multiple samples\n",
                    "       are provided for merging in this project.")
      print(ggplot() +
              annotate("text", x = 4, y = 25, size=8, label = text) +
              theme_void())
    }
  }, height = 800, width = 900)

  observeEvent(input$p17plot, {
    screenshot(id="p17plot", filename = paste("17SCA_scATAC_INTEGRATED_ATAC_SAMPLES_UMAP_TOP_MARKER_IN_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
  })

  output$p18plot <- renderPlot({
    if(length(grep("ggplot", class(data()[['p18plots']][[1]]), ignore.case = T)) > 0){
      showModal(modalDialog("Plotting figure 18..", footer=NULL, easyClose = T))
      p <- data()[['p18plots']][[paste("CLUSTER_",input$p18id, sep = "")]]
      removeModal()
    }else{
      text <- paste("\n   No multiple samples\n",
                    "       are provided for merging in this project.")
      p <- ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }
    return(p)
  }, height = 700, width = 900)

  observeEvent(input$p18plot, {
    screenshot(id="p18plot", filename = paste("18SCA_scATAC_INTEGRATED_ATAC_SAMPLES_TOP_10_MARKERS_CLUSTER_",input$p18id, sep = ""), scale = 2)
  })

  output$p19plot <- renderPlot({
    if(length(grep("ggplot", class(data()[['p19plots']][[1]]), ignore.case = T)) > 0){
      showModal(modalDialog("Plotting figure 19..", footer=NULL, easyClose = T))
      p <- data()[['p19plots']][[paste("CLUSTER_",input$p192id, sep = "")]]
      removeModal()
    }else{
      text <- paste("\n   No multiple samples\n",
                    "       are provided for merging in this project.")
      p <- ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }
    return(p)
  }, height = 800, width = 900)

  observeEvent(input$p19plot, {
    screenshot(id="p19plot", filename = paste("19SCA_scATAC_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP_PEAK_GENE_ACTIVITY_CLUSTER_",input$p19id, sep = ""), scale = 2)
  })


})


