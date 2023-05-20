###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: Spatial Transcriptomics Analysis Pipeline
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last-update date: 2023-05-11
# All Rights Reserved
###########################################################################################################################

library(BiocManager)
options(repos = BiocManager::repositories())

options(shiny.maxRequestSize=30000*1024^2)

library("shiny")
library("shinyjs")
library("shinyWidgets")
library("shinyscreenshot")
library("SingleR")
library("DT")
library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library("ggpubr")
library("ggthemes")
library("gridExtra")
library("reshape2")
library("org.Hs.eg.db")
library("xlsx")
library("ggrepel")
library("scales")
library("hdf5r")
library("SingleCellExperiment")
library("limma")

hpca.se <- readRDS("DB/hpca.se.RDS")
hs <- org.Hs.eg.db
hgnc.table <- readRDS("DB/hgnc.table.RDS")
p_val_adj <- 0.1

source("DB/SCA_Spatial_RShiny_Functions_V1.0.0.R")
color_conditions <- color_ini()

example1 <- "DB/SCA_Spatial_Example_From_10X.zip"
example2 <- "DB/SCA_Spatial_Metadata_Example.csv"

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all_metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "Percent_Mito")
p_thresh <- 0.1

ctime <- format(Sys.time(), format = "%Y%m%d%H%M%S", tz = "Europe/Stockholm")

shinyServer(function(input, output, session) {
  values <- reactiveValues(proceed = 1)
  output$keepAlive <- renderText({
    req(input$count)
    paste("keep alive ", input$count)
  })
  observeEvent(input$submit & input$termscheck,{
    
    if(input$termscheck == FALSE & input$submit){
      showModal(modalDialog("Please agree to our terms and conditions before you proceed", footer=NULL, easyClose = T))
      values$proceed <- 0
    }
    
    if(input$termscheck == TRUE & input$submit){
      if(values$proceed == 1){
        updateTabsetPanel(session, "nav",selected = "qcmetrics")
      }else{
        values$proceed <- 1
      }
    }
  })
  
  output$downloadExample1 <- downloadHandler(
    filename = "SCA_Spatial_Example_From_10X.zip",
    content = function(file) {
      file.copy(example1, file)
    },
    contentType = "application/zip"
  )
  
  output$downloadExample2 <- downloadHandler(
    filename = "SCA_Spatial_Metadata_Example.csv",
    content = function(file) {
      file.copy(example2, file)
    },
    contentType = "text/csv"
  )
  
  data <- reactive({
    req(input$files)
    showModal(modalDialog("Processing..depending on the number of samples as well as file sizes, time of processing may varies. Please wait patiently for results to be delivered.", footer=NULL))
    inFile <- input$files
    pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "SPATIAL", isDir = T)
    for(i in 1:nrow(pheno_data)){
      if(length(grep("\\.csv$|\\.txt$", pheno_data$FILE[i], ignore.case = T)) == 0){
        if(length(grep("\\.h5$", pheno_data$FILE[i], ignore.case = T)) == 0){
          pheno_data$FILE[i] <- paste(pheno_data$FILE[i],".h5", sep = "")
        }
      }
    }
    
    print(input$files)
    print(input$files$datapath)
    print(input$project_name)
    print("Project name:")
    project_name <- ifelse(input$project_name == "","My_Project",gsub("\\s+|\\(|\\)|-|\\/|\\?","",input$project_name))
    print(project_name)
    cdir <- paste(getwd(),"/",project_name,ctime, "/",sep = "")
    dir.create(cdir)
    sample_files <- unzip(inFile$datapath, exdir = cdir)
    sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
    sample_files <- sample_files[grep("\\.h5$", sample_files, ignore.case = T)]
    print("files:")
    print(sample_files)
    data <- NULL
    sample_features <- NULL
    data_current <- NULL
    data_markers <- NULL
    plotx <- NULL
    results <- NULL
    annot_names <- NULL
    scrna_data <- NULL
    removeModal()
    
    for(i in 1:length(sample_files)){
      current_sample <- gsub(".*\\/(.*\\.h5)$","\\1",sample_files[i], ignore.case = T)
      cpheno <- pheno_data[which(toupper(pheno_data$FILE) == toupper(current_sample)),]
      cname <- cpheno$SAMPLE_ID
      print(current_sample)
      showModal(modalDialog(paste("Running sample: ", current_sample, "..", sep = ""), footer=NULL))
      current <- Load10X_Spatial(data.dir = gsub("(.*\\/).*","\\1",sample_files[i], ignore.case = T), filename = current_sample)
      current <- SCTransform(current, assay = "Spatial", verbose = FALSE)
      current@project.name <- project_name
      current$orig.ident <- cname
      names(current@images) <- cname
      colnames(current@meta.data)[grep("nCount_Spatial", colnames(current@meta.data), ignore.case = T)] <- "Spatial_Counts"
      colnames(current@meta.data)[grep("nFeature_Spatial", colnames(current@meta.data), ignore.case = T)] <- "Feature_Counts"
      levels(current@active.ident) <- cname
      current <- add_names(current, cname, cname)
      current$CELL_ID <- row.names(current@meta.data)
      print(current)
      
      if(i==1){
        data <- current
        DefaultAssay(data) <- "SCT"
        VariableFeatures(data) <- VariableFeatures(data)
        sample_features[[i]] <- VariableFeatures(data)
        names(sample_features)[i] <- cname
      }else{
        data <- merge(data, current)
        DefaultAssay(data) <- "SCT"
        VariableFeatures(data) <- c(VariableFeatures(data),VariableFeatures(current))
        sample_features[[i]] <- VariableFeatures(current)
        names(sample_features)[i] <- cname
      }
      data_current[[i]] <- current
      names(data_current)[i] <- cname
      annot_names <- c(annot_names, cname)
    removeModal()
    }
    
    showModal(modalDialog(paste("Running dimension reduction and clustering..", sep = ""), footer=NULL))
    VariableFeatures(data) <- unlist(lapply(data_current, VariableFeatures))
    DefaultAssay(data) <- "SCT"
    data <- RunPCA(data, assay = "SCT", verbose = FALSE)
    data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
    data <- FindClusters(data, verbose = FALSE)
    data <- RunUMAP(data, reduction = "pca", dims = 1:30)
    # data <- RunTSNE(data, dims = 1:30, check_duplicates = FALSE)
    
    for(i in 1:length(annot_names)){
      current <- subset(data, subset = orig.ident == annot_names[i])
      current@images <- current@images[names(current@images) == annot_names[i]]
      current_clusters <- as.numeric(as.character(unique(current@meta.data$seurat_clusters)))
      Idents(current) <- "seurat_clusters"
      DefaultAssay(current) <- "Spatial"
      current_de_markers <- NULL
      current_de_markers <- FindAllMarkers(current, min.pct = 0.25, logfc.threshold = 0.25)
      current_de_markers <- data.frame(SAMPLE = annot_names[i], current_de_markers)
      current_de_markers <- current_de_markers[which(current_de_markers$p_val_adj < 0.1),]
      current_de_markers <- current_de_markers[order(current_de_markers$avg_log2FC, decreasing = T),]
      data_markers[[i]] <- current_de_markers
      names(data_markers)[i] <- annot_names[i]
      
      Idents(current) <- "seurat_clusters"
      DefaultAssay(current) <- "Spatial"
      clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(current)),
                         clusters =  current$seurat_clusters,
                         ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.main)
      current$Cell_Type <- clu_ann$labels[match(current$seurat_clusters,row.names(clu_ann))]
      current@meta.data[which(is.na(current$Cell_Type)),"Cell_Type"] <- "Unidentifiable"
      data@meta.data[which(data$orig.ident == annot_names[i]),"Cell_Type"] <- current$Cell_Type[match(row.names(data@meta.data[which(data$orig.ident == annot_names[i]),]), row.names(current@meta.data))]
      Idents(current) <- current$Cell_Type
      data_current[[which(names(data_current) %in% annot_names[i])]]$seurat_clusters <- current@meta.data[match(current$CELL_ID, data_current[[which(names(data_current) %in% annot_names[i])]]$CELL_ID),"seurat_clusters"] 
    }
    removeModal()
    # input$scrnaseq$datapath <- "/Users/luadmpan/Sites/OUTPUT/pans_test_20210822031415/Targeted_Visium_Human_Glioblastoma_Pan_Cancer_filtered_feature_bc_matrix.h5"
    # plotx <- NULL
    predictions_assay <- NULL
    scrna_files <- input$scrnaseq$datapath[grep("\\.h5$", input$scrnaseq$datapath, ignore.case = T)]
    rds_files <- input$scrnaseq$datapath[grep("\\.rds$", input$scrnaseq$datapath, ignore.case = T)]
    
    if((length(scrna_files) > 0 | length(rds_files) > 0) & length(sample_files) == 1){
      if(length(scrna_files) > 0){
        showModal(modalDialog(paste("Running processing for scRNA-Seq data.."), footer=NULL))
        scrna_data <- Read10X_h5(filename = scrna_files)
        if((length(scrna_data) == 1) | (length(scrna_data) > 10)){
          scrna_data <- CreateSeuratObject(counts = scrna_data, project = project_name, min.cells = 3, min.features = 200)
        }else{
          scrna_data <- CreateSeuratObject(counts = scrna_data$`Gene Expression`, project = prefix_name, min.cells = 3, min.features = 200)
        }
        scrna_data[["percent.mt"]] <- PercentageFeatureSet(scrna_data, pattern = "^MT-")
        scrna_data <- subset(scrna_data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
        
      }else if(length(rds_files) > 0){
        scrna_data <- readRDS(rds_files)
      }
      
      scrna_data <- SCTransform(scrna_data, ncells = 3000, verbose = FALSE)
      scrna_data <- RunPCA(scrna_data, verbose = FALSE)
      scrna_data <- RunUMAP(scrna_data, reduction = "pca", dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30))
      scrna_data <- FindNeighbors(scrna_data, reduction = "pca", dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30))
      scrna_data <- FindClusters(scrna_data, verbose = FALSE)
      
      Idents(scrna_data) <- "seurat_clusters"
      DefaultAssay(scrna_data) <- "RNA"
      clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(scrna_data)),
                         clusters =  scrna_data$seurat_clusters,
                         ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.main)
      
      scrna_data$Cell_Type <- clu_ann$labels[match(scrna_data$seurat_clusters,row.names(clu_ann))]
      scrna_data@meta.data[which(is.na(scrna_data$Cell_Type)),"Cell_Type"] <- "Unidentifiable"
      Idents(scrna_data) <- "seurat_clusters"
      DefaultAssay(scrna_data) <- "SCT"
      DefaultAssay(data) <- "Spatial"
      Idents(data) <- "seurat_clusters"
      
      data_anchors <- FindTransferAnchors(reference = scrna_data, query = data, normalization.method = "SCT")
      predictions_assay <- TransferData(anchorset = data_anchors, refdata = scrna_data$Cell_Type, 
                                        prediction.assay = TRUE, 
                                        weight.reduction = data[["pca"]], dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30))
      # data[["predictions"]] <- predictions_assay
      
      # DefaultAssay(data) <- "Spatial"
      # Idents(data) <- "Cell_Type"
      # plotx <- gen10x_plotx(scrna_data, include_meta = T)
      removeModal()
    }
    
    showModal(modalDialog("Get ready for plotting..", footer=NULL))
    results$data <- data
    if(!is.null(scrna_data)){
      results$scrna_data <- scrna_data
    }else{
      results$scrna_data <- ""
    }
    if(!is.null(predictions_assay)){
      results$predictions_assay <- predictions_assay
    }else{
      results$predictions_assay <- ""
    }
    results$data_current <- data_current
    results$sample_features <- sample_features
    results$data_markers <- data_markers
    results$cluster_colors <- gen_colors(color_conditions$tenx, length(unique(data$seurat_clusters)))
    names(results$cluster_colors) <- sort(unique(data$seurat_clusters), decreasing = F)
    results$ct_colors <- gen_colors(color_conditions$colorful, length(unique(data$Cell_Type)))
    names(results$ct_colors) <- sort(unique(data$Cell_Type))
    results$sample_colors <- gen_colors(color_conditions$bright, length(unique(data$orig.ident)))
    names(results$sample_colors) <- sort(unique(data$orig.ident), decreasing = F)
    results$project_name <- project_name
    system(paste("rm -r ",cdir, sep = ""))
    
    annot_names <- names(results$data_current)
    annot_names <- sort(annot_names)
    allclusters <- sort(unique(results$data$seurat_clusters), decreasing = F)
    # results$data_medianexpr <- data_medianexpr
    removeModal()
    
    updateSelectInput(session, inputId = 'p1id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p2id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p5id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p6id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p7id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p8id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p9id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = "p9cluster", label = "Choose a cluster to view", choices = allclusters, selected = allclusters[1])
    updateSelectInput(session, inputId = "p10id", label = "Choose a sample to display", choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = "p11id", label = "Choose a sample to display", choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = "p12id", label = "Choose a sample to display", choices = annot_names, selected = annot_names[1])
    
    # p9clusters <- sort(as.numeric(as.character(unique(results$data_current[[input$p9id]]$seurat_clusters), decreasing = F)))
    # updateSliderInput(session, inputId = 'p9cluster', label = 'Choose a cluster to view', value = p9clusters,
    #                   min = min(p9clusters), max = max(p9clusters), step = 1)
    # p9cluster = sort(unique(data()[['data_markers']][[input$p9id]]$cluster, decreasing = F))
    # print(p9cluster)
    
    # print(input$p9id)
    # updateSelectInput(session, "p9cluster", label = "Choose a cluster to view", choices = p9cluster, selected = p9cluster[1])
    print("Completed!")
    return(results)
  })
  
  output$p1plot <- renderPlot({ #renderPlotly
    # print("Data:")
    # print(data())
    plotx <- data()[['data_current']][[which(names(data()[['data_current']]) == input$p1id)]]
    plotx$orig.ident <- ifelse(nchar(plotx$orig.ident) > 15, substr(plotx$orig.ident, 1, 15), plotx$orig.ident)
    # plotx <- gen10x_plotx(plotx, selected = c("pca","umap"), include_meta = T)
    if(nrow(plotx)>0){
      showModal(modalDialog("Plotting figure 1..", footer=NULL))
      p1 <- VlnPlot(plotx, group.by = "orig.ident", features = "Spatial_Counts", pt.size = 0.1, cols = color_conditions$general) + NoLegend() + xlab("SAMPLE_NAME")
      p2 <- SpatialFeaturePlot(plotx, features = "Spatial_Counts") + theme(legend.position = "right")
      # p3 <- SpatialFeaturePlot(current, features = "Feature_Counts") + theme(legend.position = "right")
      removeModal()
      print((p1|p2)+plot_annotation(title = gsub("\\."," - ",input$p1id), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    }
  }, height = 500, width = 900)
  
  observeEvent(input$p1plot, {
    screenshot(id="p1plot", filename = paste("1SCA_VIOLIN_IMAGE_FEATURES_",input$p1id, sep = ""), scale = 2)
  })
  
  output$p2plot <- renderPlot({
    plotx <- data()[['data_current']][[which(names(data()[['data_current']]) == input$p2id)]]
    cfeatures <- data()[['sample_features']][[which(names(data()[['sample_features']]) == input$p2id)]]
    # print(cfeatures)
    DefaultAssay(plotx) <- "SCT"
    if(nrow(plotx)>0){
      showModal(modalDialog("Plotting figure 2..", footer=NULL))
      p1 <- SpatialFeaturePlot(plotx, features = head(cfeatures, 5), ncol = 5)
      # p2 <- SpatialFeaturePlot(plotx, features = head(cfeatures, 5), pt.size.factor = 1, ncol = 5)
      p2 <- SpatialFeaturePlot(plotx, features = head(cfeatures, 5), alpha = c(0.1, 1), ncol = 5)
      removeModal()
      print(p1/p2)
    }
  }, height = 500, width = 900)
  
  observeEvent(input$p2plot, {
    screenshot(id="p2plot", filename = paste("2SCA_SPATIALLY_VARIABLE_FEATURES_",input$p2id, sep = ""), scale = 2)
  })
  
  output$p3plot <- renderPlot({
    showModal(modalDialog("Plotting figure 3..", footer=NULL))
    plotx <- data()[['data']]
    plotx <- data.frame(UMAP_1 = plotx@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = plotx@reductions$umap@cell.embeddings[,"UMAP_2"],
                        CLUSTERS = plotx@meta.data$seurat_clusters,
                        SLIDES = plotx$orig.ident)
    # plotx$orig.ident <- gsub("(.*)\\..*$","\\1",plotx$orig.ident)
    
    p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTERS", point_size = 2, label_size = 6, numeric = T, annot = T,
                       plot_title = paste(data()[['project_name']], ": UMAP CLUSTERS", sep = ""), col = data()[['cluster_colors']])
    p2 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "SLIDES", point_size = 2,
                       plot_title = paste(data()[['project_name']], ": SLIDES", sep = ""), col = data()[['sample_colors']], annot = FALSE)
    p2 <- p2+theme(legend.text = element_text(size = 10))
    removeModal()
    print(p1|p2)
  }, height = 513, width = 1400)
  
  observeEvent(input$p3plot, {
    screenshot(id="p3plot", filename = paste("3SCA_UMAP_ALL_SAMPLES_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p4plot <- renderPlot({
    showModal(modalDialog("Plotting figure 4..", footer=NULL))
    plotx <- data()[['data']]
    plotx <- data.frame(UMAP_1 = plotx@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = plotx@reductions$umap@cell.embeddings[,"UMAP_2"],
                        CLUSTERS = plotx@meta.data$seurat_clusters,
                        SLIDES = plotx$orig.ident)
    # plotx$orig.ident <- gsub("(.*)\\..*$","\\1",plotx$orig.ident)
    
    p1 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", isfacet = T, title = paste(data()[['project_name']], ": UMAP CLUSTERS", sep = ""),color_by = "CLUSTERS", group_by = "SLIDES", xlabel = "UMAP_1", ylabel = "UMAP_2", strip_size = 10, col = data()[['cluster_colors']], facet_col = 2)
    p2 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", isfacet = T, title = paste(data()[['project_name']], ": SLIDES", sep = ""),color_by = "SLIDES", group_by = "SLIDES", xlabel = "UMAP_1", ylabel = "UMAP_2", strip_size = 10, col = data()[['sample_colors']], facet_col = 2)
    p2 <- p2+theme(legend.text = element_text(size = 10))
    removeModal()
    grid.arrange(p1,p2,ncol=2,widths=c(0.45, 0.55))
  }, height = 466, width = 1400)
  
  observeEvent(input$p4plot, {
    screenshot(id="p4plot", filename = paste("4SCA_UMAP_ALL_SLIDES_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p5plot <- renderPlot({
    showModal(modalDialog("Plotting figure 5..", footer=NULL))
    p <- SpatialDimPlot(data()[['data']],ncol = 1,pt.size.factor = 1.6,
                        images = input$p5id, cols = data()[['cluster_colors']])+
      ggtitle(input$p5id) +
      labs(fill= "CLUSTERS")+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.size = unit(0.5, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(fill=guide_legend(override.aes = list(size = 5)))
    removeModal()
    print(p)
  }, height = 700, width = 900)
  
  observeEvent(input$p5plot, {
    screenshot(id="p5plot", filename = paste("5SCA_UMAP_SLIDE_CLUSTERS_",input$p5id, sep = ""), scale = 2)
  })
  
  output$p6plot <- renderPlot({
    showModal(modalDialog("Plotting figure 6..", footer=NULL))
    current <- data()[['data_current']][[which(names(data()[['data_current']]) == input$p6id)]]
    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "SCT"
    col_left <- as.character(unique(current$seurat_clusters))[match(levels(data()[['data']]@meta.data$seurat_clusters), as.character(unique(current$seurat_clusters)))]
    current@images <- current@images[names(current@images) == input$p6id]
    removeModal()
    print(SpatialDimPlot(current, cells.highlight = CellsByIdentities(object = current,
                                                                      idents = as.numeric(as.character(col_left[!is.na(col_left)]))), facet.highlight = TRUE, ncol = 6))
  }, height = 600, width = 900)
  
  observeEvent(input$p6plot, {
    screenshot(id="p6plot", filename = paste("6SCA_SLIDE_IMAGE_BY_EACH_CLUSTER_",input$p6id, sep = ""), scale = 2)
  })
  
  output$p7plot <- renderPlot({
    showModal(modalDialog("Plotting figure 7..", footer=NULL))
    plotx <- subset(data()[['data']], subset = orig.ident == input$p7id)
    Idents(plotx) <- "seurat_clusters"
    plotx <- data.frame(UMAP_1 = plotx@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = plotx@reductions$umap@cell.embeddings[,"UMAP_2"],
                        CLUSTERS = plotx$seurat_clusters)
    removeModal()
    print(plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTERS", plot_title = paste(data()[['project_name']], ": UMAP CLUSTERS", sep = ""), annot = T, label_size = 6, point_size = 2, numeric = T, col = data()[['cluster_colors']]))
  }, height = 600, width = 900)
  
  observeEvent(input$p7plot, {
    screenshot(id="p7plot", filename = paste("7SCA_UMAP_SAMPLE_CLUSTER_",input$p7id, sep = ""), scale = 2)
  })
  
  output$p8plot <- renderPlot({
    showModal(modalDialog("Plotting figure 8..", footer=NULL))
    plotx <- subset(data()[['data']], subset = orig.ident == input$p8id)
    Idents(plotx) <- "seurat_clusters"
    plotx <- data.frame(Cluster = plotx$seurat_clusters, 
                        X = plotx@images[[input$p8id]]@coordinates$imagerow,
                        Y = plotx@images[[input$p8id]]@coordinates$imagecol)
    
    col_left <- as.character(unique(plotx$Cluster))[match(levels(data()[['data']]@meta.data$seurat_clusters), as.character(unique(plotx$Cluster)))]
    
    p1 <- plot_slice(plotx, pt_size = 3, plot_title = "SLICE CLUSTER VIEW", col = data()[['cluster_colors']][which(!is.na(col_left))], is.facet = FALSE, annot = FALSE)
    p1 <- p1+theme(legend.title = element_text(size = 15),
                   legend.text = element_text(size = 15),
                   legend.key.size = unit(0.2, "cm"))
    p2 <- plot_slice(plotx, pt_size = 1.5, plot_title = "SLICE CLUSTER VIEW", col = data()[['cluster_colors']][which(!is.na(col_left))], is.facet = TRUE, annot = FALSE, strip_size = 15)
    removeModal()
    print((p1|p2)+plot_annotation(title = input$p8id, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  }, height = 450, width = 900)
  
  observeEvent(input$p8plot, {
    screenshot(id="p8plot", filename = paste("8SCA_UMAP_SAMPLE_IMAGE_CLUSTER_",input$p8id, sep = ""), scale = 2)
  })
  
  output$p9plot <- renderPlot({
    showModal(modalDialog("Plotting figure 9..", footer=NULL))
    current <- subset(data()[['data']], subset = orig.ident == input$p9id)
    current_de <- data()[['data_markers']][[input$p9id]]
    current_de <- current_de[which(as.character(current_de$cluster) == as.character(input$p9cluster)),]
    if(nrow(current_de) > 0){
      current_de <- current_de[1:ifelse(nrow(current_de) < 5, nrow(current_de), 5),]
      p <- SpatialFeaturePlot(object = current, features = current_de$gene, images = input$p9id, ncol = ifelse(nrow(current_de) < 5, nrow(current_de), 5))
    }else{
      print("current_de")
      print(current_de)
      text <- paste("\n   Sample does not belong to this cluster")
      p <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
    }
    removeModal()
    return(p)
  }, height = 300, width = 900)
  
  observeEvent(input$p9plot, {
    screenshot(id="p9plot", filename = paste("9SCA_SPATIALLY_TOP_FEATURES_",input$p9id, sep = ""), scale = 2)
  })
  
  options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
  output$p10table <- renderDT(data()[['data_markers']][[input$p10id]],
                              filter = "top",
                              style="bootstrap",
                              rownames = F,
                              options = list(pageLength = 10))
  
  output$p11plot <- renderPlot({
    showModal(modalDialog("Plotting figure 11..", footer=NULL))
    current <- subset(data()[['data']], subset = orig.ident == input$p11id)
    current_de_markers <- data()[['data_markers']][[input$p11id]]
    current_de_markers <- current_de_markers[current_de_markers$p_val_adj < 0.1,]
    current_de_markers <- current_de_markers[order(current_de_markers$avg_log2FC, decreasing = T),]
    current_de_markers <- split(current_de_markers, current_de_markers$cluster)
    top_markers <- NULL
    for(j in 1:length(current_de_markers)){
      if(nrow(current_de_markers[[j]]) > 0){
        top_markers <- rbind(top_markers,current_de_markers[[j]][1,])
      }
    }
    
    DefaultAssay(current) <- "Spatial"
    p <- VlnPlot(current, group.by = "seurat_clusters", features = unique(top_markers$gene),
                 pt.size = 0, ncol = 2, cols = data()[['cluster_colors']], log = T)
    p <- p+plot_annotation(title = paste("TOP GENE IN EACH CLUSTER\nLog(Average Expression)", sep = ""), theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
    removeModal()
    return(p)
  }, height = 1100, width = 900)
  
  observeEvent(input$p11plot, {
    screenshot(id="p11plot", filename = paste("11SCA_VIOLINPLOT_TOP_MARKERS_",input$p11id, sep = ""), scale = 2)
  })
  
  output$p12plot <- renderPlot({
    showModal(modalDialog("Plotting figure 12..", footer=NULL))
    current <- subset(data()[['data']], subset = orig.ident == input$p12id)
    Idents(current) <- current$Cell_Type
    plotx <- data.frame(UMAP_1 = current@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = current@reductions$umap@cell.embeddings[,"UMAP_2"],
                        CELL_TYPE = current$Cell_Type)
    p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = "SPATIAL UMAP - CELL TYPE",
                       col = data()[['ct_colors']], annot = TRUE, legend_position = "right", point_size = 2)
    p2 <- SpatialDimPlot(current,ncol = 1,pt.size.factor = 1.6,images = input$p12id,
                         cols = data()[['ct_colors']][which(names(data()[['ct_colors']]) %in% unique(current$Cell_Type))])+
      ggtitle("SPATIAL SLIDE - CELL TYPE") +
      labs(fill= "CELL TYPES")+
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1, "cm"))+
      guides(fill=guide_legend(override.aes = list(size = 5)))
    removeModal()
    return((p1/p2)+plot_annotation(title = paste(data()[['project_name']],": ",input$p12id,sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  }, height = 1000, width = 900)
  
  observeEvent(input$p12plot, {
    screenshot(id="p12plot", filename = paste("12SCA_UMAP_AUTO_CELLTYPE_ANNOTATIONS_",input$p12id, sep = ""), scale = 2)
  })
  
  output$p13plot <- renderPlot({
    if(class(data()[['scrna_data']])[1] == "Seurat"){
      showModal(modalDialog("Plotting figure 13..", footer=NULL))
      plotx <- data()[['scrna_data']]
      Idents(plotx) <- "seurat_clusters"
      plotx <- data.frame(UMAP_1 = plotx@reductions$umap@cell.embeddings[,"UMAP_1"],
                          UMAP_2 = plotx@reductions$umap@cell.embeddings[,"UMAP_2"],
                          CLUSTERS = plotx$seurat_clusters)
      p <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTERS", plot_title = paste("scRNASEQ - UMAP CLUSTERS: ", data()[['project_name']], sep = ""), annot = T, label_size = 8, point_size = 2, numeric = T, col = gen_colors(color_conditions$monet, length(unique(plotx$CLUSTERS))))
      removeModal()
    }else{
      text <- paste("\n   No scRNA-Seq data\n",
                    "       is provided in this project.")
      p <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
    }
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p13plot, {
    screenshot(id="p13plot", filename = paste("13SCA_UMAP_scRNASEQ_SAMPLE_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p14plot <- renderPlot({
    if(class(data()[['scrna_data']])[1] == "Seurat"){
      showModal(modalDialog("Plotting figure 14..", footer=NULL))
      plotx <- data()[['scrna_data']]
      Idents(plotx) <- "Cell_Type"
      plotx <- data.frame(UMAP_1 = plotx@reductions$umap@cell.embeddings[,"UMAP_1"],
                          UMAP_2 = plotx@reductions$umap@cell.embeddings[,"UMAP_2"],
                          CELL_TYPE = plotx$Cell_Type)
      p <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = paste("scRNASEQ - UMAP CELL TYPES: ", data()[['project_name']], sep = ""), annot = T, label_size = 8, point_size = 2, numeric = F, col = gen_colors(color_conditions$vanngogh, length(unique(plotx$CELL_TYPE))))
      removeModal()
    }else{
      text <- paste("\n   No scRNA-Seq data\n",
                    "       is provided in this project.")
      p <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
    }
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p14plot, {
    screenshot(id="p14plot", filename = paste("14SCA_UMAP_scRNASEQ_CELL_TYPE_ANNOTATION_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p15plot <- renderPlot({
    if(class(data()[['scrna_data']])[1] == "Seurat"){
      showModal(modalDialog("Plotting figure 15..", footer=NULL))
      plotx <- data()[['data']]
      plotx[['predictions']] <- data()[['predictions_assay']]
      DefaultAssay(plotx) <- "predictions"
      Idents(plotx) <- "Cell_Type"
      p <- VlnPlot(plotx, group.by = "seurat_clusters", features = row.names(plotx)[grep("^max$", row.names(plotx), ignore.case = T, invert = T)],pt.size = 0, ncol = 2, cols = data()[['cluster_colors']])
      p <- p+plot_annotation(title = paste("CELL TYPE PREDICTION SCORES FOR EACH SPATIAL CLUSTER", sep = ""), theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))&xlab("CLUSTERS")
      removeModal()
    }else{
      text <- paste("\n   No scRNA-Seq data\n",
                    "       is provided in this project.")
      p <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
    }
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p15plot, {
    screenshot(id="p15plot", filename = paste("15SCA_VIOLIN_PLOT_scRNASEQ_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p16plot <- renderPlot({
    if(class(data()[['scrna_data']])[1] == "Seurat"){
      showModal(modalDialog("Plotting figure 16..", footer=NULL))
      plotx <- data()[['data']]
      plotx[['predictions']] <- data()[['predictions_assay']]
      DefaultAssay(plotx) <- "predictions"
      Idents(plotx) <- "Cell_Type"
      p <- SpatialFeaturePlot(plotx, features = row.names(plotx)[grep("^max$", row.names(plotx), ignore.case = T, invert = T)], pt.size.factor = 1.6, ncol = 1, crop = TRUE)
      p <- p+plot_annotation(title = paste("CELL TYPE PREDICTION SCORE: ", data()[['project_name']], sep = ""), theme = theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)))
      removeModal()
    }else{
      text <- paste("\n   No scRNA-Seq data\n",
                    "       is provided in this project.")
      p <- ggplot() + 
        annotate("text", x = 4, y = 25, size=8, label = text) + 
        theme_void()
    }
    return(p)
  }, width = 1400)
  
  observeEvent(input$p16plot, {
    screenshot(id="p16plot", filename = paste("16SCA_UMAP_scRNASEQ_PREDICTED_CELL_TYPE_ANNOTATION_",data()[['project_name']], sep = ""), scale = 2)
  })
})
