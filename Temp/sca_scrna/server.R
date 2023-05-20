###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: 10X scRNA Sequencing Analysis Pipeline
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last-update date: 2023-05-11
# All Rights Reserved
###########################################################################################################################
library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=30000*1024^2)

library("Seurat")
library("monocle3")
library("dplyr")
library("patchwork")
library("SeuratWrappers")
library("reshape2")
library("reticulate")
library("plotly")
library("viridis")
library("cowplot")
library("plot3D")
library("DDRTree")
library("pheatmap")
library("scater")
library("ggbeeswarm")
library("ggthemes")
library("gplots")
library("gridExtra")
library("ggpubr")
library("plotly")
library("lattice")
library("ReactomePA")
library("HGNChelper")
library("DOSE")
library("enrichplot")
library("harmony")
library("ggrepel")
library("akmedoids")
library("grid")
library("multtest")
library("metap")
library("akmedoids")
library("htmlwidgets")
library("celldex")
library("ggnewscale")
library("clusterProfiler")
library("ggupset")
library("europepmc")
library("dendextend")
library("org.Hs.eg.db")
library("scRNAseq")
library("SingleR")
library("DT")
library("R.utils")
library("hdf5r")
library("ggridges")
library("scales")
library("xlsx")

hs <- org.Hs.eg.db
hgnc.table <- readRDS("DB/hgnc.table.RDS")

Sys.setenv("DISPLAY"=":0")

source("DB/SCA_scRNASEQ_RShiny_Functions_V1.0.0.R")
color_conditions <- color_ini()

hpca.se <- readRDS("DB/HumanPrimaryCellAtlasData.RDS")

example1 <- "DB/SCA_scRNASEQ_Example_From_10X.zip"
example2 <- "DB/SCA_scRNASEQ_Metadata_Example.csv"

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

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
    
    output$downloadExample1 <- downloadHandler(
        filename = "SCA_scRNASEQ_Example_From_10X.zip",
        content = function(file) {
            file.copy(example1, file)
        },
        contentType = "application/zip"
    )
    
    output$downloadExample2 <- downloadHandler(
        filename = "SCA_scRNASEQ_Metadata_Example.csv",
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
        pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "scRNASEQ", isDir = T)
        for(i in 1:nrow(pheno_data)){
            if(length(grep("\\.csv$|\\.txt$", pheno_data$FILE[i], ignore.case = T)) == 0){
                if(length(grep("\\.h5$", pheno_data$FILE[i], ignore.case = T)) == 0){
                    pheno_data$FILE[i] <- paste(pheno_data$FILE[i],".h5", sep = "")
                }
            }
        }
        pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, sep = "_")
        sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SID)))
        names(sample_colors) <- unique(pheno_data$SID)
        
        print(input$files)
        print(input$files$datapath)
        project_name <- ifelse(input$project_name == "","My_Project",gsub("\\s+|\\(|\\)|-|\\/|\\?","",input$project_name))
        print(project_name)
        print("==========================================")
        print("==========================================")
        print(paste("Initiating new project: ", project_name, ctime, sep = ""))
        print("==========================================")
        print("==========================================")
        
        cdir <- paste(getwd(),"/",project_name,ctime, "/",sep = "")
        dir.create(cdir)
        
        sample_files <- unzip(input$files$datapath, exdir = cdir)
        sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
        
        removeModal()
        
        #######################################################################################################################################
        showModal(modalDialog("Initializing..", footer=NULL, easyClose = F))
        data <- NULL
        data_current <- NULL
        annot_names <- NULL
        results <- NULL
        removeModal()
        
        for(j in 1:nrow(pheno_data)){
            showModal(modalDialog(paste("Running ", pheno_data[j,"FILE"], "..", sep = ""), footer=NULL, easyClose = F))
            annot_names <- c(annot_names, pheno_data[j,"SID"])
            data_current[[j]] <- NULL
            data_current[[j]] <- sample_files[grep(pheno_data[j,"FILE"], sample_files, ignore.case = T)]
            if(length(grep("\\.h5$", data_current[[j]], ignore.case = T)) > 0){
                data_current[[j]] <- Read10X_h5(data_current[[j]])
            }else{
                data_current[[j]] <- readfile(data_current[[j]])
            }
            if((length(data_current[[j]]) == 1) | (length(data_current[[j]]) > 10)){
                data_current[[j]] <- CreateSeuratObject(counts = data_current[[j]], project = annot_names[j], min.cells = 3, min.features = 200)
            }else{
                data_current[[j]] <- CreateSeuratObject(counts = data_current[[j]]$`Gene Expression`, project = annot_names[j], min.cells = 3, min.features = 200)
            }
            
            data_current[[j]][["Percent_Mito"]] <- PercentageFeatureSet(data_current[[j]], pattern = "^MT-")
            if(max(data_current[[j]][["Percent_Mito"]]) == 0){
                plot_mito <- FALSE
            }else{
                plot_mito <- TRUE
            }
            
            names(data_current)[j] <- annot_names[j]
            data_current[[j]] <- add_names(data_current[[j]], sample_name = annot_names[j], current_ident = "orig.ident")
            data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, pheno_data[j,which(colnames(pheno_data) != "SAMPLE_ID")])
            plotx <- data_current[[j]]@meta.data
            colnames(plotx)[grep("orig.ident", colnames(plotx), ignore.case = T)] <- "SAMPLE_ID"
            
            p1 <- NULL
            p2 <- NULL
            p <- NULL
            p1 <- own_violin(plotx, feature = "nFeature_RNA", plotx_title = "No. Genes Detected / Cell", col = color_conditions$tenx[1], title.size = 18)
            p2 <- own_violin(plotx, feature = "nCount_RNA", plotx_title = "No. Molecules Detected / Cell", col = color_conditions$tenx[2], title.size = 18)
            if(plot_mito == TRUE){
                p3 <- own_violin(plotx, feature = "Percent_Mito", plotx_title = "Mitochondria Percent / Cell", col = color_conditions$tenx[3], title.size = 18)
                p <- p1+p2+p3
            }else{p <- p1+p2}
            results[['p1plots']][[j]] <- p+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
            names(results[['p1plots']])[j] <- pheno_data[j,"SID"]
            
            p1 <- NULL
            p2 <- NULL
            p <- NULL
            p2 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                              title_name = annot_names[j], col = color_conditions$tenx[5],
                              xlabel = "No. Molecules Detected/Cell", ylabel = "No. Genes Detected/Cell")
            if(plot_mito == FALSE){
                p <- p2 }else{
                    p1 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "Percent_Mito",
                                      title_name = annot_names[j], col = color_conditions$tenx[4],
                                      xlabel = "No. Molecules Detected/Cell", ylabel = "Mitochondria Percent/Cell")
                    p <- p1+p2
                }
            
            results[['p2plots']][[j]] <- p + plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5)))
            names(results[['p2plots']])[j] <- pheno_data[j,"SID"]
            
            data_current[[j]] <- NormalizeData(data_current[[j]], verbose = TRUE)
            # data_current[[j]]@assays$RNA@data@x[is.na(data_current[[j]]@assays$RNA@data@x)] <- 0
            data_current[[j]] <- FindVariableFeatures(data_current[[j]],selection.method = "vst")
            data_current[[j]] <- ScaleData(data_current[[j]], features = row.names(data_current)[j])
            
            if(nrow(data_current[[j]]) > 2000){
                orig_gene_names <- NULL
                scale_orig_gene_names <- NULL
                if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
                    orig_gene_names <- row.names(data_current)[j]
                    scale_orig_gene_names <- row.names(data_current[[j]]@assays$RNA@scale.data)
                    row.names(data_current[[j]]@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@counts), ignore.case = T)
                    row.names(data_current[[j]]@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@data), ignore.case = T)
                    row.names(data_current[[j]]@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@scale.data), ignore.case = T)
                }
                
                data_current[[j]] <- CellCycleScoring(data_current[[j]], g2m.features=g2m.genes, s.features=s.genes, set.ident = TRUE)
                data_current[[j]]@meta.data$Phase <- factor(data_current[[j]]@meta.data$Phase, levels = c("G1","S","G2M"))
                data_current[[j]] <- RunPCA(data_current[[j]], features = c(s.genes, g2m.genes))
                data_current[[j]]@reductions$pca_selected <- data_current[[j]]@reductions$pca
                
                if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
                    row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
                    row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
                    row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
                }
                
                data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
                
                # G2/M and S Phase Markers: Tirosh et al, 2015
                p1 <- own_2d_scatter(data_current[[j]], "pca", "Phase", "PCA Based on Variable Features")
                p2 <- own_2d_scatter(data_current[[j]], "pca_selected", "Phase", "PCA Based on G2/M and S Phase Markers")
                results[['p3plots']][[j]] <- p1 / p2 + plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p3plots']])[j] <- pheno_data[j,"SID"]
                
                data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca_selected", dims = 1:ifelse(length(data_current[[j]]@reductions$pca_selected) < 30, length(data_current[[j]]@reductions$pca_selected), 30))
                data_current[[j]]@reductions$umap_selected <- data_current[[j]]@reductions$umap
                data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
                
                p1 <- own_2d_scatter(data_current[[j]], "umap", "Phase", "UMAP Based on Variable Features")
                p2 <- own_2d_scatter(data_current[[j]], "umap_selected", "Phase", "UMAP Based on G2/M and S Phase Markers")
                results[['p4plots']][[j]] <- p1 / p2 + plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p4plots']])[j] <- pheno_data[j,"SID"]
                
            }else{
                data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
                data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
            }
            
            metrics <- colnames(data_current[[j]]@meta.data)[grep("nCount_RNA|nFeature_RNA|S.Score|G2M.Score|Percent_Mito", colnames(data_current[[j]]@meta.data), ignore.case = T)]
            
            results[['p5plots']][[j]] <- FeaturePlot(data_current[[j]],
                                                     reduction = "umap",
                                                     features = metrics,
                                                     pt.size = 0.4,
                                                     order = TRUE,
                                                     min.cutoff = 'q10',
                                                     label = FALSE, cols = c("green","blue"))+
                plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""),
                                theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
            names(results[['p5plots']])[j] <- pheno_data[j,"SID"]
            
            ########################################################################################################################
            # Filtering
            data_current[[j]] <- subset(data_current[[j]], subset = nFeature_RNA > 200 & nFeature_RNA <= 25000 & Percent_Mito < 5)
            if(ncol(data_current[[j]]) > 200){
                
                data_current[[j]] <- NormalizeData(data_current[[j]], verbose = TRUE)
                # data_current[[j]]@assays$RNA@data@x[is.na(data_current[[j]]@assays$RNA@data@x)] <- 0
                data_current[[j]] <- FindVariableFeatures(data_current[[j]],selection.method = "vst")
                data_current[[j]] <- ScaleData(data_current[[j]], features = rownames(data_current)[j])
                
                if(input$ccregression != "NO_REGRESSION" & (nrow(data_current[[j]][['RNA']]) > 2000)){
                    
                    orig_gene_names <- NULL
                    scale_orig_gene_names <- NULL
                    
                    if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
                        orig_gene_names <- row.names(data_current)[j]
                        scale_orig_gene_names <- row.names(data_current[[j]]@assays$RNA@scale.data)
                        row.names(data_current[[j]]@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@counts), ignore.case = T)
                        row.names(data_current[[j]]@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@data), ignore.case = T)
                        row.names(data_current[[j]]@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@scale.data), ignore.case = T)
                    }
                    
                    data_current[[j]] <- CellCycleScoring(data_current[[j]], g2m.features=g2m.genes, s.features=s.genes, set.ident = TRUE)
                    
                    data_current[[j]]@meta.data$Phase <- factor(data_current[[j]]@meta.data$Phase, levels = c("G1","S","G2M"))
                    
                    if((input$ccregression == "TWO_PHASES") | (input$ccregression == "YES")){
                        cc_method <- c("S.Score", "G2M.Score")
                        cc_type <- "Phase Scores"
                        
                    }else if (input$ccregression == "PHASE_DIFFERENCE"){
                        data_current[[j]]$CC.Difference <- data_current[[j]]$S.Score - data_current[[j]]$G2M.Score
                        cc_method <- "CC.Difference"
                        cc_type <- "Difference"
                    }
                    
                    data_current[[j]] <- ScaleData(data_current[[j]], vars.to.regress = cc_method, features = rownames(data_current)[j])
                    data_current[[j]] <- RunPCA(data_current[[j]], features = c(s.genes, g2m.genes))
                    data_current[[j]]@reductions$pca_selected <- data_current[[j]]@reductions$pca
                    
                    if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
                        row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
                        row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
                        row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
                    }
                    
                    data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
                    
                    p1 <- DimPlot(data_current[[j]], reduction = "pca", split.by = "Phase", cols = color_conditions$bright) +
                        ggtitle(paste("PCA Based on Variable Features", sep = ""))+
                        theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
                    p2 <- DimPlot(data_current[[j]], reduction = "pca_selected", split.by = "Phase", cols = color_conditions$bright) +
                        ggtitle(paste("PCA Based on G2/M and S Phase Markers", sep = ""))+
                        theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
                    
                    data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca_selected", dims = 1:ifelse(length(data_current[[j]]@reductions$pca_selected) < 30, length(data_current[[j]]@reductions$pca_selected), 30))
                    data_current[[j]]@reductions$umap_selected <- data_current[[j]]@reductions$umap
                    data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
                    
                    metrics <- colnames(data_current[[j]]@meta.data)[grep("nCount_RNA|nFeature_RNA|S.Score|G2M.Score|Percent_Mito", colnames(data_current[[j]]@meta.data), ignore.case = T)]
                    
                    p1 <- NULL
                    p2 <- NULL
                    p1 <- DimPlot(data_current[[j]], reduction = "pca",split.by = "Phase", cols = color_conditions$tenx) +
                        ggtitle(paste("PCA Based on Variable Features", sep = "")) +
                        theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
                    p2 <- DimPlot(data_current[[j]], reduction = "pca_selected", split.by = "Phase", cols = color_conditions$tenx) +
                        ggtitle(paste("PCA Based on G2/M and S Phase Markers", sep = "")) +
                        theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
                    
                    # G2/M and S Phase Markers: Tirosh et al, 2015
                    results[['p6plots']][[j]] <- p1 / p2 + plot_annotation(title = paste("Post Cell Cycle Regression: ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                    names(results[['p6plots']])[j] <- pheno_data[j,"SID"]
                    
                    p1 <- NULL
                    p2 <- NULL
                    p1 <- DimPlot(data_current[[j]], reduction = "umap", cols = color_conditions$tenx, split.by = "Phase") +
                        ggtitle(paste("UMAP Based on Variable Features", sep = "")) + theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
                    p2 <- DimPlot(data_current[[j]], reduction = "umap_selected", cols = color_conditions$tenx, split.by = "Phase") +
                        ggtitle(paste("UMAP Based on G2/M and S Phase Markers", sep = "")) + theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
                    results[['p7plots']][[j]] <- p1 / p2 + plot_annotation(title = paste("Post Cell Cycle Regression: ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                    names(results[['p7plots']])[j] <- pheno_data[j,"SID"]
                    
                    results[['p8plots']][[j]] <- FeaturePlot(data_current[[j]],
                                                             reduction = "umap",
                                                             features = metrics,
                                                             pt.size = 0.4,
                                                             order = TRUE,
                                                             min.cutoff = 'q10',
                                                             label = FALSE, cols = c("green","blue"))+
                        plot_annotation(title = paste("UMAP Post Cell Cycle Regression - ",annot_names[j], sep = ""),
                                        theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5)))
                    names(results[['p8plots']])[j] <- pheno_data[j,"SID"]
                    
                }else{
                    DefaultAssay(data_current[[j]]) <- 'RNA'
                    data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
                    data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
                    
                    metrics <- colnames(data_current[[j]]@meta.data)[grep("nCount_RNA|nFeature_RNA|S.Score|G2M.Score|Percent_Mito", colnames(data_current[[j]]@meta.data), ignore.case = T)]
                    
                    results[['p8plots']][[j]] <- FeaturePlot(data_current[[j]],
                                                             reduction = "umap",
                                                             features = metrics,
                                                             pt.size = 0.4,
                                                             order = TRUE,
                                                             min.cutoff = 'q10',
                                                             label = FALSE, cols = c("green","blue"))+
                        plot_annotation(title = paste("UMAP Features - ",annot_names[j], sep = ""),
                                        theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5)))
                    names(results[['p8plots']])[j] <- pheno_data[j,"SID"]
                    
                }
                
                colnames(data_current[[j]]@assays$RNA@counts) <- gsub("_|-|\\s+|\\t","",colnames(data_current[[j]]@assays$RNA@counts))
                colnames(data_current[[j]]@assays$RNA@data) <- gsub("_|-|\\s+|\\t","",colnames(data_current[[j]]@assays$RNA@data))
                colnames(data_current[[j]]@assays$RNA@scale.data) <- gsub("_|-|\\s+|\\t","",colnames(data_current[[j]]@assays$RNA@scale.data))
                row.names(data_current[[j]]@meta.data) <- gsub("_|-|\\s+|\\t","",row.names(data_current[[j]]@meta.data))
                
                # if(length(h5_files) > 0){
                #   SCTRANS <- TRUE
                #   current <- suppressWarnings(SCTransform(data_current[[j]], verbose = FALSE))
                # }else{
                #   SCTRANS <- FALSE
                # }
                
                plotx <- data_current[[j]]@meta.data
                colnames(plotx)[grep("orig.ident", colnames(plotx), ignore.case = T)] <- "SAMPLE_ID"
                
                p1 <- NULL
                p2 <- NULL
                p <- NULL
                p1 <- own_violin(plotx, feature = "nFeature_RNA", plotx_title = "No. Genes Detected / Cell", col = color_conditions$tenx[1], title.size = 15)
                p2 <- own_violin(plotx, feature = "nCount_RNA", plotx_title = "No. Molecules Detected / Cell", col = color_conditions$tenx[2], title.size = 15)
                if(plot_mito == TRUE){
                    p3 <- own_violin(plotx, feature = "Percent_Mito", plotx_title = "Mitochondria Percent / Cell", col = color_conditions$tenx[3], title.size = 15)
                    p <- p1+p2+p3
                }else{
                    p <- p1+p2
                }
                results[['p9plots']][[j]] <- p+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p9plots']])[j] <- pheno_data[j,"SID"]
                
                p2 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                                  title_name = annot_names[j], col = color_conditions$tenx[5],
                                  xlabel = "No. Molecules Detected/Cell", ylabel = "No. Genes Detected/Cell")
                if(plot_mito == TRUE){
                    p1 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "Percent_Mito",
                                      title_name = annot_names[j], col = color_conditions$tenx[4],
                                      xlabel = "No. Molecules Detected/Cell", ylabel = "Mitochondria Percent/Cell")
                    p <- p1+p2
                }else{
                    p <- p2
                }
                
                results[['p9plots']][[j]] <- p + plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
                names(results[['p9plots']])[j] <- pheno_data[j,"SID"]
                
                n <- 10
                p1 <- NULL
                p2 <- NULL
                p1 <- VariableFeaturePlot(data_current[[j]], cols = color_conditions$bright[1:2])
                p2 <- LabelPoints(plot = p1, points = VariableFeatures(data_current[[j]])[1:n], repel = TRUE, xnudge = 0, ynudge = 0)
                results[['p11plots']][[j]] <- p2+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
                names(results[['p11plots']])[j] <- pheno_data[j,"SID"]
                
                results[['p12plots']][[j]] <- VizDimLoadings(data_current[[j]], dims = 1:2, reduction = "pca",
                                                             col = color_conditions$manycolors[1])+
                    plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
                names(results[['p12plots']])[j] <- pheno_data[j,"SID"]
                
                Idents(data_current[[j]]) <- "orig.ident"
                data_current[[j]] <- RunPCA(data_current[[j]])
                results[['p13plots']][[j]] <- DimPlot(data_current[[j]], reduction = "pca",cols = color_conditions$tenx) + theme(legend.position = "none")+
                    plot_annotation(title = paste(annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p13plots']])[j] <- pheno_data[j,"SID"]
                
                results[['p14plots']][[j]] <- DimHeatmap(data_current[[j]], dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE, assays = "RNA")+
                    plot_annotation(title = paste(annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p14plots']])[j] <- pheno_data[j,"SID"]
                
                data_current[[j]] <- FindNeighbors(data_current[[j]], dims = 1:ifelse(ncol(data_current[[j]]) > 30, 30, ncol(data_current[[j]])))
                data_current[[j]] <- FindClusters(data_current[[j]], resolution = 0.8)
                data_current[[j]] <- RunTSNE(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30), check_duplicates = FALSE)
                data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
                
                plotx <- gen10x_plotx(data_current[[j]])
                plotx$CLUSTER <- Idents(data_current[[j]])
                plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(as.numeric(as.character(unique(plotx$CLUSTER)))))
                ccolors <- gen_colors(color_conditions$colorful,length(levels(plotx$CLUSTER)))
                names(ccolors) <- levels(plotx$CLUSTER)
                
                p1 <- NULL
                p2 <- NULL
                p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER", plot_title = "UMAP",point_size = 1,
                                   col = color_conditions$colorful, annot = TRUE, legend_position = "right", numeric = TRUE)
                p2 <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CLUSTER", plot_title = "tSNE",point_size = 1,
                                   col = color_conditions$colorful, annot = TRUE, legend_position = "right", numeric = TRUE)
                
                results[['p15plots']][[j]] <- p1+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p15plots']])[j] <- pheno_data[j,"SID"]
                
                results[['p16plots']][[j]] <- p2+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p16plots']])[j] <- pheno_data[j,"SID"]
                
                current_clusters <- sort(as.numeric(as.character(unique(plotx$CLUSTER))))
                Idents(data_current[[j]]) <- "seurat_clusters"
                current_out <- NULL
                current_out <- deanalysis(data_current[[j]], current_clusters, plot_title = annot_names[j],group=NULL,de_analysis = "findallmarkers")
                current_data_markers <- NULL
                results[['p17data']][[j]] <- current_out$current_data_markers
                names(results[['p17data']])[j] <- pheno_data[j,"SID"]
                de_type <- current_out$de_type
                de_name <- current_out$de_name
                top1 <- current_out$top1
                topn <- current_out$topn
                current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
                ######################### MEDIAN EXPRESSIONS #####################################
                current <- group_medianexpr(current_out$current_data_markers, data = data_current[[j]], group = "seurat_clusters", cell_type = F)
                plot_median <- current$plot_median
                top_markers <- current$top_markers
                plot_median[which(is.na(plot_median))] <- 0
                results[['p18plots']][[j]] <- complex_heatmap(plot_median, col = color_conditions$BlueYellowRed)
                names(results[['p18plots']])[j] <- pheno_data[j,"SID"]
                
                #################################### TOP MARKERS + TSNE / UMAP ###################################################################
                wt <- colnames(current_out$current_data_markers)[grep("log.*FC", colnames(current_out$current_data_markers), ignore.case = T)]
                current_clusters <- sort(as.numeric(as.character(unique(top_markers$cluster))), decreasing = F)
                logFC_list <- NULL
                p1 <- NULL
                P2 <- NULL
                P3 <- NULL
                plotx <- gen10x_plotx(data_current[[j]], groups = NULL)
                plotx$CLUSTER <- data_current[[j]]$seurat_clusters
                ccolors <- gen_colors(color_conditions$colorful,length(unique(plotx$CLUSTER)))
                names(ccolors) <- sort(unique(plotx$CLUSTER))
                top_markers <- data.frame(top_markers)
                top_markers$cluster <- factor(top_markers$cluster, levels = sort(unique(top_markers$cluster)))

                p1 <- ggplot(data=top_markers, aes_string(x="gene", y=wt, fill="cluster")) +
                  geom_bar(stat="identity", position=position_dodge())+
                    scale_fill_manual(values = ccolors)+
                    theme_classic() +
                    coord_flip()+
                    ylab(wt)+
                    facet_wrap(~cluster, scales = "free_y") +
                    theme(axis.text.y=element_text(size = 10),
                          strip.text.x = element_text(size = 15, face = "bold"),
                          legend.title = element_text(size =20, face = "bold"),
                          legend.text = element_text(size = 15))
                
                p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CLUSTER",
                                   plot_title = annot_names[j],col = ccolors,numeric = T,
                                   annot = T, legend_position = "right", point_size = 1)
                p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CLUSTER",
                                   plot_title = annot_names[j],col = ccolors,numeric = T,
                                   annot = T, legend_position = "right", point_size = 1)
                
                results[['p19plots']][[j]] <- p1
                names(results[['p19plots']])[j] <- pheno_data[j,"SID"]
                results[['p192plots']][[j]] <- p2
                names(results[['p192plots']])[j] <- pheno_data[j,"SID"]
                results[['p193plots']][[j]] <- p3
                names(results[['p193plots']])[j] <- pheno_data[j,"SID"]
                
                # grid.arrange(p1, p2, ncol = 2)
                # grid.arrange(p1, p3, ncol = 2)
                
                # top_1_markers <- current_out$current_data_markers %>% group_by(cluster) %>% top_n(n = 1, eval(parse(text = wt)))
                
                top_1_markers <- split(current_out$current_data_marker, current_out$current_data_marker$cluster)
                top_1_markers <- lapply(top_1_markers, function(x){
                    x <- x[order(x[,wt], decreasing = T),]
                    x <- x[1,]
                })
                top_1_markers <- do.call(rbind.data.frame, top_1_markers)
                
                results[['p20plots']][[j]] <- RidgePlot(data_current[[j]], features = unique(unlist(top_1_markers$gene)), ncol = 4,
                                                        cols = ccolors)+
                    plot_annotation(title = paste("TOP1: ", annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p20plots']])[j] <- pheno_data[j,"SID"]
                
                results[['p21plots']][[j]] <- current_out$featureplot +plot_layout(ncol = 4) +
                    plot_annotation(title = paste("TOP1: ",de_name,annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                names(results[['p21plots']])[j] <- pheno_data[j,"SID"]
                
                for(k in 1:length(current_clusters)){
                    current <- topn[which(topn$cluster == current_clusters[k]),]
                    current <- current[order(current$p_val_adj, decreasing = F),]
                    results[['p22plots']][[length(results[['p22plots']])+1]] <- VlnPlot(data_current[[j]], features = current$gene, pt.size = 0,
                                                                                        ncol = 4, cols = ccolors)&
                        xlab("CLUSTERS")&
                        plot_annotation(title = paste("TOP",n," IN CLUSTER ", current_clusters[k],": ",current_out$de_name, annot_names[j],  sep = ""),
                                        theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
                    names(results[['p22plots']])[length(results[['p22plots']])] <- paste("CLUSTER_",current_clusters[k],"|",pheno_data[j,"SID"], sep = "")
                }
                
                results[['p23plots']][[j]] <- DoHeatmap(data_current[[j]], features = topn$gene,
                                                        group.colors = ccolors, size = 8) +
                    ggtitle(annot_names[j])+
                    NoLegend()+theme(axis.text.x = element_blank(),
                                     axis.text.y = element_text(size = 10),
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                     legend.title = element_blank(),
                                     legend.text = element_text(size = 15),
                                     plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
                names(results[['p23plots']])[j] <- pheno_data[j,"SID"]
                
                Idents(data_current[[j]]) <- "seurat_clusters"
                DefaultAssay(data_current[[j]]) <- "RNA"
                
                orig_gene_names <- NULL
                scale_orig_gene_names <- NULL
                if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
                    orig_gene_names <- row.names(data_current)[j]
                    scale_orig_gene_names <- row.names(data_current[[j]]@assays$RNA@scale.data)
                    row.names(data_current[[j]]@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@counts), ignore.case = T)
                    row.names(data_current[[j]]@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@data), ignore.case = T)
                    row.names(data_current[[j]]@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@scale.data), ignore.case = T)
                }
                
                clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data_current[[j]])),
                                   clusters =data_current[[j]]$seurat_clusters,
                                   ref = hpca.se, assay.type.test=1,
                                   labels = hpca.se$label.main)
                data_current[[j]]$CELL_TYPE <- clu_ann$labels[match(data_current[[j]]$seurat_clusters,row.names(clu_ann))]
                data_current[[j]]@meta.data[which(is.na(data_current[[j]]$CELL_TYPE)),"CELL_TYPE"] <- "Unidentifiable"
                if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
                    row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
                    row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
                    row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
                }
                
                Idents(data_current[[j]]) <- data_current[[j]]$CELL_TYPE
                plotx <- gen10x_plotx(data_current[[j]])
                plotx$CELL_TYPE <- data_current[[j]]$CELL_TYPE
                
                cct_colors <- gen_colors(color_conditions$tenx,length(unique(plotx$CELL_TYPE)))
                names(cct_colors) <- sort(unique(plotx$CELL_TYPE))
                
                results[['p24plots']][[j]] <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = annot_names[j],
                                                           col = cct_colors, annot = TRUE, legend_position = "right", point_size = 1)
                names(results[['p24plots']])[j] <- pheno_data[j,"SID"]
                
                results[['p25plots']][[j]] <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CELL_TYPE", plot_title = annot_names[j],
                                                           col = cct_colors, annot = TRUE, legend_position = "right", point_size = 1)
                names(results[['p25plots']])[j] <- pheno_data[j,"SID"]
                
                ################ Pathway Analysis #########################################################################
                filtered_markers <- current_out$current_data_markers[grep("AC[0-9]+\\.[0-9]+|AL[0-9]+\\.[0-9]+",current_out$current_data_markers$gene, ignore.case = T, invert = T),]
                topn <- filtered_markers[which(abs(filtered_markers$avg_log2FC) > 0.25),]
                topn <- topn[order(topn$avg_log2FC, decreasing = T),]
                
                if(length(grep("ENS.*-.*", topn$gene)) > length(topn$gene)/2){
                    topn$ENSEMBL_ID <- gsub("(ENS.*?[0-9]+)[-|_|\\s+|\\.].*","\\1",topn$gene)
                    mapped_id <- select(hs, keys = topn$ENSEMBL_ID, columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
                    topn$ENTREZ_ID <- mapped_id[match(topn$ENSEMBL_ID, mapped_id$ENSEMBL),"ENTREZID"]
                }else{
                    mapped_id <- select(hs, keys = topn$gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
                    topn$ENTREZ_ID <- mapped_id[match(topn$gene, mapped_id$SYMBOL),"ENTREZID"]
                }
                
                if(length(which(is.na(topn$ENTREZ_ID))) > 0 & !(length(grep("ENS.*-.*", topn$gene)) > length(topn$gene)/2)){
                    current <- topn[which(is.na(topn$ENTREZ_ID)),]
                    current$alternate_gene <- hgnc.table[match(current$gene, hgnc.table$Symbol),"Approved.Symbol"]
                    if(!all(current$gene == current$alternate_gene)){
                      current_id <- select(hs, keys = current$alternate_gene[!is.na(current$alternate_gene)], columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
                      current$ENTREZ_ID <- current_id[match(current$alternate_gene, current_id$SYMBOL),"ENTREZID"]
                      topn[which(is.na(topn$ENTREZ_ID)),"ENTREZ_ID"] <- current[match(as.character(unlist(topn[which(is.na(topn$ENTREZ_ID)),"gene"])), current$gene),"ENTREZ_ID"]
                    }
                }
                topn <- topn[!is.na(topn$ENTREZ_ID),]
                
                p_threshold <- 0.01
                current <- split(topn, topn$cluster)
                pathway_EA_result <- NULL
                pathway_EA_result <- lapply(current, function(x){
                    x <- enrichDGN(gene=unique(as.numeric(as.character(x$ENTREZ_ID))),pvalueCutoff=0.05, qvalueCutoff = 0.05,readable=T)
                })
                
                for(i in 1:length(pathway_EA_result)){
                    if(!is.null(pathway_EA_result[[i]])){
                        pathway_EA_result[[i]]@result <- pathway_EA_result[[i]]@result[which(pathway_EA_result[[i]]@result$p.adjust < p_threshold),]
                    }
                }
                
                results[['pathway_EA_result']][[j]] <- pathway_EA_result
                names(results[['pathway_EA_result']])[j] <- pheno_data[j,"SID"]
                
                summary_pathways <- NULL
                n <- 15
                for(i in 1:length(pathway_EA_result)){
                    current <- pathway_EA_result[[i]]
                    if(!is.null(current)){
                        if(nrow(current@result) > 0){
                            summary_pathways <- rbind(summary_pathways, data.frame(CLUSTER = names(pathway_EA_result)[i],current))
                            results[['p26plots']][[length(results[['p26plots']])+1]] <- barplot(current, showCategory=n, orderBy = "x")+
                                ggtitle(paste(annot_names[j],"\nEnriched Terms for ORA: Cluster ", names(pathway_EA_result)[i],sep = ""))
                            names(results[['p26plots']])[length(results[['p26plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                        }
                    }
                }
                
                results[['p27data']][[j]] <- summary_pathways
                names(results[['p27data']])[j] <- pheno_data[j,"SID"]
                
                for(i in 1:length(pathway_EA_result)){
                    current <- pathway_EA_result[[i]]
                    if(!is.null(current)){
                        if(nrow(current@result) > 0){
                            results[['p28plots']][[length(results[['p28plots']])+1]] <- dotplot(current, showCategory=n, orderBy = "x")+ggtitle(paste(annot_names[j], "\nEnriched Terms for ORA: Cluster ", names(pathway_EA_result)[i],sep = ""))
                            names(results[['p28plots']])[length(results[['p28plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                        }
                    }
                }
                
                topn_cluster <- split(topn, topn$cluster)
                p1 <- NULL
                p2 <- NULL
                p3 <- NULL
                
                for(i in 1:length(pathway_EA_result)){
                    current <- pathway_EA_result[[i]]
                    if(!is.null(current)){
                        if(nrow(current@result) > 0){
                            current_gene <- topn_cluster[[i]]
                            if(toupper(wt) == toupper("avg_log2FC")){
                                gene_list <- current_gene$avg_log2FC
                            }else{
                                gene_list <- current_gene$avg_logFC
                            }
                            names(gene_list) <- current_gene$ENTREZ_ID
                            gene_list <- sort(gene_list, decreasing = TRUE)
                            gene_list <- gene_list[!duplicated(names(gene_list))]
                            
                            p1[[i]] <- cnetplot(current, foldChange=gene_list) +
                                ggtitle(paste(annot_names[j], "\nA: Gene-Concept Network, ",wt,": Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                theme(plot.title = element_text(size = 6, face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
                            # p2[[i]] <- cnetplot(current, categorySize="pvalue", foldChange=gene_list) +
                            #     ggtitle(paste(annot_names[j], "\nB: Gene-Concept Network, ",wt," & P-Value: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                            #     theme(plot.title = element_text(size = 6, face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
                            p2[[i]] <- cnetplot(current, foldChange=gene_list, circular = TRUE, colorEdge = TRUE) +
                                ggtitle(paste(annot_names[j], "\nC: Gene-Concept Network, ",wt,": Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                theme(plot.title = element_text(size = 6, face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
                            
                            results[['p29plots']][[length(results[['p29plots']])+1]] <- cowplot::plot_grid(p1[[i]], p2[[i]], ncol=2, rel_widths=c(1.2, 1.2, 1.2))
                            names(results[['p29plots']])[length(results[['p29plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                        }
                    }
                }
                
                p1 <- NULL
                p2 <- NULL
                
                for(i in 1:length(pathway_EA_result)){
                    current <- pathway_EA_result[[i]]
                    if(!is.null(current)){
                        if(nrow(current@result) > 0){
                            p1[[i]] <- cnetplot(current, node_label="category") +
                                ggtitle(paste(annot_names[j], "\nA: Network Terms: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
                            p2[[i]] <- cnetplot(current, node_label="gene") +
                                ggtitle(paste(annot_names[j], "\nA: Network Genes: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
                            
                            results[['p30plots']][[length(results[['p30plots']])+1]] <- cowplot::plot_grid(p1[[i]], p2[[i]], ncol=2, rel_widths=c(1.2, 1.2, 1.2))
                            names(results[['p30plots']])[length(results[['p30plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                        }
                    }
                }
                
                for(i in 1:length(pathway_EA_result)){
                    current <- pathway_EA_result[[i]]
                    if(!is.null(current)){
                        if(nrow(current@result) > 0){
                            current_gene <- topn_cluster[[i]]
                            if(toupper(wt) == toupper("avg_log2FC")){
                                gene_list <- current_gene$avg_log2FC
                            }else{
                                gene_list <- current_gene$avg_logFC
                            }
                            names(gene_list) <- current_gene$ENTREZ_ID
                            gene_list <- sort(gene_list, decreasing = TRUE)
                            gene_list <- gene_list[!duplicated(names(gene_list))]
                            # currentR <- setReadable(current, 'org.Hs.eg.db', 'ENTREZID')
                            results[['p31plots']][[length(results[['p31plots']])+1]] <- heatplot(current, foldChange = gene_list) +
                                ggtitle(paste(annot_names[j], "\nHeatmap of Enrichment Terms: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                theme(plot.title = element_text(face = "bold"))
                            names(results[['p31plots']])[length(results[['p31plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                        }
                    }
                }
                
                for(i in 1:length(pathway_EA_result)){
                    current <- pathway_EA_result[[i]]
                    if(!is.null(current)){
                        if(nrow(current@result) > 0){
                            currentR <- pairwise_termsim(current)
                            if(nrow(currentR@result) > 10){
                                results[['p32plots']][[length(results[['p32plots']])+1]] <- emapplot(currentR, cex_line = 0.1) +
                                    ggtitle(paste(annot_names[j], "\nEnrichment Map: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                    scale_colour_gradient_tableau(palette = "Red-Gold") +
                                    theme(plot.title = element_text(face = "bold"))
                                names(results[['p32plots']])[length(results[['p32plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                            }
                        }
                    }
                }
                
                current <- NULL
                for(i in 1:length(topn_cluster)){
                    if(!is.null(topn_cluster[[i]])){
                        current[[i]] <- topn_cluster[[i]]$ENTREZ_ID
                    }
                }
                
                names(current) <- names(topn_cluster)
                topn_cluster <- topn_cluster[!unlist(lapply(topn_cluster,is.null))]
                current <- current[!unlist(lapply(current,is.null))]
                plotx <- compareCluster(current, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)
                currentR <- pairwise_termsim(plotx)
                
                results[['p33plots']][[j]] <- emapplot(currentR, cex_line = 0.1, pie="count", pie_scale=1.5, layout="kk") +
                    ggtitle(paste(annot_names[j], "\nEnrichment Map with Proportions", sep = "")) +
                    scale_fill_tableau(palette = "Tableau 20") +
                    theme(plot.title = element_text(face = "bold"))
                names(results[['p33plots']])[j] <- pheno_data[j,"SID"]
                
                for(i in 1:length(pathway_EA_result)){
                    if(!is.null(pathway_EA_result[[i]])){
                        if(nrow(pathway_EA_result[[i]]@result) > 0){
                            results[['p34plots']][[length(results[['p34plots']])+1]] <- upsetplot(pathway_EA_result[[i]]) +
                                ggtitle(paste(annot_names[j], "\nUpset Plot: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
                            names(results[['p34plots']])[length(results[['p34plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                        }
                    }
                }
                
                GSEA_result <- NULL
                for(i in 1:length(topn_cluster)){
                    if(toupper(wt) == toupper("avg_logFC")){
                        gene_list <- topn_cluster[[i]]$avg_logFC
                    }else{
                        gene_list <- topn_cluster[[i]]$avg_log2FC
                    }
                    if(nrow(topn_cluster[[i]])>10){
                        names(gene_list) <- topn_cluster[[i]]$ENTREZ_ID
                        gene_list <- sort(gene_list, decreasing = TRUE)
                        gene_list <- gene_list[!duplicated(names(gene_list))]
                        GSEA_result[[i]] <- gseNCG(gene=gene_list, pvalueCutoff = p_threshold)
                    }else{
                        GSEA_result[[i]] <- NULL
                    }
                }
                
                if_plot <- NULL
                
                for(i in 1:length(GSEA_result)){
                    if(!is.null(GSEA_result[[i]])){
                        if(nrow(GSEA_result[[i]]@result) > 0){
                            GSEA_result[[i]]@result <- GSEA_result[[i]]@result[GSEA_result[[i]]@result$p.adjust < p_threshold,]
                            if_plot[[i]] <- ifelse(nrow(GSEA_result[[i]]@result) ==0, FALSE, TRUE)
                        }else{
                            if_plot[[i]] <- FALSE
                        }
                    }else{
                        if_plot[[i]] <- FALSE
                    }
                }
                
                if(!all(if_plot == FALSE)){
                    
                    for(i in 1:length(GSEA_result)){
                        current <- GSEA_result[[i]]
                        if(!is.null(current)){
                            if(nrow(current@result[which(current@result$ID != "-"),]) > 0){
                                results[['p35plots']][[length(results[['p35plots']])+1]] <- gseaplot2(current, geneSetID = which(current@result$ID != "-"), pvalue_table = TRUE,
                                                                                                      ES_geom = "line") # , title = current$Description[which(current$Description != "-")], 
                                #+
                                    # ggtitle(paste(annot_names[j], "\nGSEA Plot: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                                    # theme(plot.title = element_text(face = "bold"))
                                names(results[['p35plots']])[length(results[['p35plots']])] <- paste("CLUSTER_",names(pathway_EA_result)[i],"|",pheno_data[j,"SID"],sep = "")
                            }
                        }
                    }
                }
                
                ##################### CELL TYPES ##########################################################################
                
                if(length(unique(data_current[[j]]$CELL_TYPE)) > 1){
                    current <- group_medianexpr(current_out$current_data_markers, data = data_current[[j]], group = "CELL_TYPE", cell_type = T)
                    plot_median_cell_type <- current$plot_median
                    top_markers_cell_type <- current$top_markers
                    n <- length(unique(top_markers_cell_type$gene))
                    plot_median_cell_type[which(is.na(plot_median_cell_type))] <- 0
                    results[['p36plots']][[j]] <- complex_heatmap(plot_median_cell_type, col = color_conditions$BlueYellowRed)
                    names(results[['p36plots']])[j] <- pheno_data[j,"SID"]
                    
                    wt <- colnames(current_out$current_data_markers)[grep("log.*FC", colnames(current_out$current_data_markers), ignore.case = T)]
                    
                    plotx <- gen10x_plotx(data_current[[j]], groups = NULL)
                    plotx$CELL_TYPE <- data_current[[j]]$CELL_TYPE
                    
                    # cell_types <- sort(as.character(unique(plotx$CELL_TYPE)))
                    # colors <- gen_colors(color_conditions$colorful,length(cell_types))
                    cct_colors <- gen_colors(color_conditions$tenx,length(unique(plotx$CELL_TYPE)))
                    names(cct_colors) <- sort(unique(plotx$CELL_TYPE))
                    top_markers_cell_type <- data.frame(top_markers_cell_type)
                    top_markers_cell_type$CELL_TYPE <- factor(top_markers_cell_type$CELL_TYPE, levels = sort(unique(top_markers_cell_type$CELL_TYPE)))
                    p1 <- NULL
                    p2 <- NULL
                    p3 <- NULL
                    p1 <- ggplot(data=top_markers_cell_type, aes_string(x="gene", y=wt, fill="CELL_TYPE")) +
                      geom_bar(stat="identity", position=position_dodge())+
                        scale_fill_manual(values = cct_colors)+
                        theme_classic() +
                        coord_flip()+
                        ylab(wt)+
                        facet_wrap(~CELL_TYPE, scales = "free_y") +
                        theme(axis.text.y=element_text(size = 10),
                              strip.text.x = element_text(size = 15, face = "bold"),
                              legend.title = element_text(size =20, face = "bold"),
                              legend.text = element_text(size = 15))
                    
                    p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CELL_TYPE",
                                       plot_title = annot_names[j],col = cct_colors,numeric = F,
                                       annot = T, legend_position = "right", point_size = 1)
                    p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CELL_TYPE",
                                       plot_title = annot_names[j],col = cct_colors,numeric = F,
                                       annot = T, legend_position = "right", point_size = 1)
                    
                    results[['p37plots']][[j]] <- p1
                    names(results[['p37plots']])[j] <- pheno_data[j,"SID"]
                    results[['p372plots']][[j]] <- p2
                    names(results[['p372plots']])[j] <- pheno_data[j,"SID"]
                    results[['p373plots']][[j]] <- p3
                    names(results[['p373plots']])[j] <- pheno_data[j,"SID"]
                    
                    # grid.arrange(p1, p2, ncol = 2)
                    # grid.arrange(p1, p3, ncol = 2)
                    
                    current_out$current_data_markers$CELL_TYPE <- data_current[[j]]@meta.data[match(current_out$current_data_markers$cluster, data_current[[j]]$seurat_clusters),"CELL_TYPE"]
                    # top_1_markers <- current_out$current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = 1, eval(parse(text = wt)))
                    top_1_markers <- split(current_out$current_data_marker, current_out$current_data_marker$CELL_TYPE)
                    top_1_markers <- lapply(top_1_markers, function(x){
                        x <- x[order(x[,wt], decreasing = T),]
                        x <- x[1,]
                    })
                    top_1_markers <- do.call(rbind.data.frame, top_1_markers)
                    Idents(data_current[[j]]) <- "CELL_TYPE"
                    results[['p38plots']][[j]] <- RidgePlot(data_current[[j]], features = unique(unlist(top_1_markers$gene)), ncol = 2,
                                                            cols = cct_colors)+
                        plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
                    names(results[['p38plots']])[j] <- pheno_data[j,"SID"]
                    
                    Idents(data_current[[j]]) <- "orig.ident"
                    DefaultAssay(data_current[[j]]) <- 'RNA'
                    rm(p)
                }
            }
        }
        ########## Integration Steps #####################################################################
        if(nrow(pheno_data) == 1){
            
            data <- data_current[[1]]
            DefaultAssay(data) <- "RNA"
            n <- 10
            top_features <- data.frame(TOP_PCA_POS_NEG_GENES = paste(capture.output(print(data[["pca"]], dims = 1:5, nfeatures = n))))
            
        }else{
            data_current <- lapply(data_current, function(x){
                # x <- ScaleData(x, verbose=F, features = data_features, vars.to.regress = c("nCount_RNA", "Percent_Mito"))
                # x <- RunPCA(x, npcs = 30, verbose = FALSE, features = data_features)
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                # x <- RunPCA(x, npcs = 30, verbose = FALSE, features = VariableFeatures(x))
            })
            integrative_features <- SelectIntegrationFeatures(object.list = data_current)
            data_anchors <- FindIntegrationAnchors(object.list = data_current,
                                                   reduction = "rpca", anchor.features = integrative_features)
            data <- IntegrateData(anchorset = data_anchors)
            rm(data_anchors)
            DefaultAssay(data) <- "integrated"
            data <- ScaleData(data, verbose = FALSE)
            
            DefaultAssay(data) <- "RNA"
            data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
            data <- ScaleData(data, verbose = FALSE)
            data <- RunPCA(data, verbose = FALSE)
            data <- RunUMAP(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30))
            data <- RunTSNE(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30), check_duplicates = FALSE)
            data_dim <- data.frame(gen10x_plotx(data), DATA_TYPE = "BEFORE_INTEGRATION", SAMPLE_ID = data$orig.ident)
            
            if((toupper(input$integration_method) == "SEURAT") | (toupper(input$integration_method) == "NULL")){
                # integration_method <- "SEURAT"
                integration_name <- "SEURAT_INTEGRATED"
                DefaultAssay(data) <- "integrated"
                reduction_method <- "pca"
                
                data <- RunPCA(data, verbose = FALSE)
                data <- RunUMAP(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30))
                data <- RunTSNE(data, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30), check_duplicates = FALSE)
                current <- cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident))
            }else if(toupper(input$integration_method) == "HARMONY"){
                integration_name <- "HARMONY_INTEGRATED"
                DefaultAssay(data) <- "RNA"
                reduction_method <- "harmony"
                data <- RunHarmony(data, group.by.vars = "BATCH")
                data <- RunUMAP(data, reduction = "harmony", dims = 1:ifelse(length(data@reductions$harmony) < 30, length(data@reductions$harmony), 20))
                data <- RunTSNE(data, reduction = "harmony", dims = 1:ifelse(length(data@reductions$harmony) < 30, length(data@reductions$harmony), 20), check_duplicates = FALSE)
                current <- cbind(data.frame(genharmony_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident))
            }
            
            data_dim <- rbind(data_dim, cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident)))
            
            results[['p39plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                         data_dim[data_dim$DATA_TYPE == integration_name,],
                                                         dim1= "UMAP_1", dim2 = "UMAP_2", group = "SAMPLE_ID",
                                                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                         maintitle = project_name, titlesize = 35, col = sample_colors)
            
            results[['p40plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                         data_dim[data_dim$DATA_TYPE == integration_name,],
                                                         dim1= "tSNE_1", dim2 = "tSNE_2", group = "SAMPLE_ID",
                                                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                         maintitle = project_name, titlesize = 35, col = sample_colors)
            
            if(input$integration_method == "SEURAT"){
                
                results[['p41plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                             data_dim[data_dim$DATA_TYPE == integration_name,],
                                                             dim1= "PC_1", dim2 = "PC_2", group = "SAMPLE_ID",
                                                             subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                             maintitle = project_name, titlesize = 35, col = sample_colors)
            }else if(input$integration_method == "HARMONY"){
                p1 <- NULL
                p2 <- NULL
                p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",], x = "PC_1", y = "PC_2", group = "SAMPLE_ID", plot_title = "BEFORE_INTEGRATION",
                                   col = sample_colors, annot = FALSE, legend_position = "bottom")
                p2 <- plot_bygroup(current, x = "HARMONY_1", y = "HARMONY_2", group = "SAMPLE_ID", plot_title = integration_name,
                                   col = sample_colors, annot = FALSE, legend_position = "bottom")
                
                results[['p41plots']] <- p1+p2+
                    plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5)))
            }
            
            data_dim$GROUP <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"GROUP"]
            data_dim$BATCH <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"BATCH"]
            
            cgroup_colors <- gen_colors(color_conditions$bright,length(unique(data_dim$GROUP)))
            names(cgroup_colors) <- sort(unique(data_dim$GROUP))
            
            results[['p42plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                         data_dim[data_dim$DATA_TYPE == integration_name,],
                                                         dim1= "UMAP_1", dim2 = "UMAP_2",group = "GROUP",
                                                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                         maintitle = project_name, titlesize = 35, col = cgroup_colors)
            
            
            results[['p43plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                         data_dim[data_dim$DATA_TYPE == integration_name,],
                                                         dim1= "tSNE_1", dim2 = "tSNE_2",group = "GROUP",
                                                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                         maintitle = project_name, titlesize = 35, col = cgroup_colors)
            
            p1 <- NULL
            p2 <- NULL
            p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",], x = "PC_1", y = "PC_2", group = "GROUP", plot_title = "BEFORE_INTEGRATION",
                               col = cgroup_colors, annot = FALSE, legend_position = "bottom", point_size = 1)
            if(input$integration_method == "SEURAT"){
                p2 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == integration_name,], x = "PC_1", y = "PC_2", group = "GROUP", plot_title = integration_name,
                                   col = cgroup_colors, annot = FALSE, legend_position = "bottom", point_size = 1)
            }else{
                p2 <- plot_bygroup(current, x = "HARMONY_1", y = "HARMONY_2", group = "GROUP", plot_title = integration_name,
                                   col = cgroup_colors, annot = FALSE, legend_position = "bottom", point_size = 1)
                
            }
            results[['p44plots']] <- p1+p2+
                plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5)))
            
            cbatch_colors <- gen_colors(color_conditions$general,length(unique(data_dim$BATCH)))
            names(cbatch_colors) <- sort(unique(data_dim$BATCH))
            
            results[['p45plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                         data_dim[data_dim$DATA_TYPE == integration_name,],
                                                         dim1= "UMAP_1", dim2 = "UMAP_2", group = "BATCH",
                                                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                         maintitle = project_name, titlesize = 35, col = cbatch_colors)
            
            results[['p46plots']] <- beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                         data_dim[data_dim$DATA_TYPE == integration_name,],
                                                         dim1= "tSNE_1", dim2 = "tSNE_2", group = "BATCH",
                                                         subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                         maintitle = project_name, titlesize = 35, col = cbatch_colors)
            
            p1 <- NULL
            p2 <- NULL
            p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",], x = "PC_1", y = "PC_2", group = "BATCH", plot_title = "BEFORE_INTEGRATION",point_size = 1,
                               col = cbatch_colors, annot = FALSE, legend_position = "bottom")
            if(input$integration_method == "SEURAT"){
                p2 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == integration_name,], x = "PC_1", y = "PC_2", group = "BATCH", plot_title = integration_name,point_size = 1,
                                   col = cbatch_colors, annot = FALSE, legend_position = "bottom")
            }else{
                p2 <- plot_bygroup(current, x = "HARMONY_1", y = "HARMONY_2", group = "BATCH", plot_title = integration_name,
                                   col = cbatch_colors, point_size = 1,annot = FALSE, legend_position = "bottom")
            }
            results[['p47plots']] <- p1+p2+
                plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5)))
            
            
            DefaultAssay(data) <- "integrated"
            
            n <- 10
            top_features <- data.frame(TOP_PCA_POS_NEG_GENES = paste(capture.output(print(data[[ifelse(input$integration_method == "SEURAT", "pca", "harmony")]], dims = 1:5, nfeatures = n))))
            
            results[['p48plots']] <- VizDimLoadings(data, dims = 1:2, reduction = ifelse(input$integration_method == "SEURAT", "pca", "harmony"),col = color_conditions$mark[1])+ plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
            
            ############## Integration Plots #################################################################
            DefaultAssay(data) <- "integrated"
            results[['p49plots']] <- DimHeatmap(data, assays = ifelse(input$integration_method == "SEURAT", "integrated", "RNA"),
                                                reduction = ifelse(input$integration_method == "SEURAT", "pca", "harmony"),
                                                dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)+
                plot_annotation(title = paste(project_name, sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
            
            if(input$integration_method == "SEURAT"){
              DefaultAssay(data) <- "integrated"
            }else{
              DefaultAssay(data) <- "RNA"
            }
            
            data <- FindNeighbors(data, reduction = reduction_method, dims = 1:ifelse(ncol(data) < 30, ncol(data), 30))
            data <- FindClusters(data, resolution = 0.8)
            
            if(input$integration_method == "HARMONY"){
                integration_cluster <- "RNA_snn_res.0.8"
            }else{
                integration_cluster <- "integrated_snn_res.0.8"
            }
            
            plotx <- gen10x_plotx(data)
            plotx$CLUSTER <- Idents(data)
            plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(as.numeric(as.character(unique(plotx$CLUSTER)))))
            plotx$GROUP <- data$GROUP
            current_clusters <- sort(as.numeric(as.character(unique(Idents(data)))),decreasing = F)
            
            groups <- sort(unique(plotx$GROUP))
            cluster_colors <- gen_colors(color_conditions$colorful,length(levels(plotx$CLUSTER)))
            names(cluster_colors) <- levels(plotx$CLUSTER)
            
            p1 <- NULL
            p2 <- NULL
            p1 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", title = project_name,isfacet = T,
                                    col=cluster_colors, color_by = "CLUSTER", group_by = "GROUP",
                                    xlabel = "UMAP_1", ylabel = "UMAP_2")
            p2 <- own_facet_scatter(plotx, "tSNE_1", "tSNE_2", title = project_name,isfacet = T,
                                    col=cluster_colors, color_by = "CLUSTER", group_by = "GROUP",
                                    xlabel = "UMAP_1", ylabel = "UMAP_2")
            
            results[['p50plots']] <- p1
            results[['p51plots']] <- p2
            
            Idents(data) <- integration_cluster
            DefaultAssay(data) <- "RNA"
            current_out <- deanalysis(data, current_clusters, plot_title = project_name,group=NULL,de_analysis = "findallmarkers", cluster_name = integration_cluster)
            results[['p52data']] <- current_out$current_data_markers
            de_type <- current_out$de_type
            de_name <- current_out$de_name
            top1 <- current_out$top1
            topn <- current_out$topn
            current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
            
            results[['p53plots']] <- current_out$featureplot+plot_layout(ncol = 4) +
                # facet_wrap(~GROUP)+
                plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5)))
            
            for(k in 1:length(current_clusters)){
                current <- topn[which(topn$cluster == current_clusters[k]),]
                current <- current[order(current$p_val_adj, decreasing = F),]
                results[['p54plots']][[length(results[['p54plots']])+1]] <- VlnPlot(data, features = unique(current$gene), pt.size = 0, split.by = "GROUP",
                                                                                    stack = T,flip = T,
                                                                                    cols = cgroup_colors)&
                    xlab("CLUSTERS")&
                    plot_annotation(title = paste(project_name, ": CLUSTER ", current_clusters[k], sep = ""),
                                    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
                names(results[['p54plots']])[length(results[['p54plots']])] <- current_clusters[k]
            }
            
            results[['p55plots']] <- DoHeatmap(data, features = topn$gene,
                                               group.colors = cluster_colors, size = 8) +
                ggtitle(project_name)+
                NoLegend()+theme(axis.text.x = element_blank(),
                                 axis.text.y = element_text(size = 10),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                 legend.title = element_blank(),
                                 legend.text = element_text(size = 15),
                                 plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
            
            DefaultAssay(data) <- "RNA"
            orig_gene_names <- NULL
            scale_orig_gene_names <- NULL
            if(length(grep("ENSG[0-9]+",row.names(data), ignore.case = T)) > nrow(data)/2){
                orig_gene_names <- row.names(data)
                scale_orig_gene_names <- row.names(data@assays$RNA@scale.data)
                row.names(data@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data@assays$RNA@counts), ignore.case = T)
                row.names(data@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data@assays$RNA@data), ignore.case = T)
                row.names(data@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data@assays$RNA@scale.data), ignore.case = T)
            }
            
            Idents(data) <- "seurat_clusters"
            DefaultAssay(data) <- "RNA"
            clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
                               clusters =  data$seurat_clusters,
                               ref = hpca.se, assay.type.test=1,
                               labels = hpca.se$label.main)
            data$CELL_TYPE <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
            data@meta.data[which(is.na(data$CELL_TYPE)),"CELL_TYPE"] <- "Unidentifiable"
            
            if(length(grep("ENSG[0-9]+",row.names(data), ignore.case = T)) > nrow(data)/2){
                row.names(data@assays$RNA@counts) <- orig_gene_names
                row.names(data@assays$RNA@data) <- orig_gene_names
                row.names(data@assays$RNA@scale.data) <- scale_orig_gene_names
            }
            Idents(data) <- data$CELL_TYPE
            plotx <- gen10x_plotx(data)
            plotx$CELL_TYPE <- data$CELL_TYPE
            plotx$GROUP <- data$GROUP
            
            ct_colors <- gen_colors(color_conditions$bright,length(unique(plotx$CELL_TYPE)))
            names(ct_colors) <- sort(unique((plotx$CELL_TYPE)))
            
            results[['p56plots']] <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = project_name,
                                                  col = ct_colors, annot = TRUE, legend_position = "right", point_size = 1)
            
            results[['p57plots']] <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CELL_TYPE", plot_title = project_name,
                                                  col = ct_colors, annot = TRUE, legend_position = "right", point_size = 1)
            
            results[['p58plots']] <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", isfacet = T,
                                                       color_by = "CELL_TYPE", xlabel = "UMAP_1", ylabel = "UMAP_2",
                                                       group_by = "GROUP", title = project_name,
                                                       col = ct_colors)
            
            #################################### MEDIAN EXPRESSION #########################################################
            
            current <- group_medianexpr(current_out$current_data_markers, data, ref_group = integration_cluster, group = "seurat_clusters", cell_type = F)
            plot_median <- current$plot_median
            top_markers <- current$top_markers
            plot_median[which(is.na(plot_median))] <- 0
            results[['p59plots']] <- complex_heatmap(plot_median, col = color_conditions$BlueYellowRed)
            
            #################################### TOP MARKERS + TSNE / UMAP ###################################################################
            current_clusters <- sort(as.numeric(as.character(unique(top_markers$cluster))), decreasing = F)
            wt <- colnames(top_markers)[grep("log.*FC", colnames(top_markers), ignore.case = T)]
            
            p1 <- NULL
            p2 <- NULL
            p3 <- NULL
            top_markers <- data.frame(top_markers)
            top_markers$cluster <- factor(top_markers$cluster, levels = current_clusters)
            p1 <- ggplot(data=top_markers, aes_string(x="gene", y=wt, fill="cluster")) +
                geom_bar(stat="identity", position=position_dodge())+
                scale_fill_manual(values = cluster_colors)+
                theme_classic() +
                coord_flip()+
                ylab(wt)+
                facet_wrap(~cluster, scales = "free_y") +
                theme(axis.text.y=element_text(size = 10),
                      strip.text.x = element_text(size = 15, face = "bold"),
                      legend.title = element_text(size =20, face = "bold"),
                      legend.text = element_text(size = 15))
            
            plotx <- gen10x_plotx(data, groups = NULL)
            plotx$CLUSTER <- data@meta.data[,integration_cluster]
            plotx$GROUP <- data$GROUP
            
            p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CLUSTER",
                               plot_title = project_name,col = cluster_colors,numeric = T,
                               annot = T, legend_position = "right", point_size = 3)
            p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CLUSTER",
                               plot_title = project_name,col = cluster_colors,numeric = T,
                               annot = T, legend_position = "right", point_size = 3)
            
            results[['p60plots']] <- p1
            results[['p602plots']] <- p2
            results[['p603plots']] <- p3
            
            # grid.arrange(p1, p2, ncol = 2)
            # grid.arrange(p1, p3, ncol = 2)
            
            current_clusters <- sort(as.numeric(as.character(unique(top_markers$cluster))), decreasing = F)
            # top_1_markers <- current_out$current_data_markers %>% group_by(cluster) %>% top_n(n = 1, eval(parse(text = wt)))
            
            top_1_markers <- split(current_out$current_data_marker, current_out$current_data_marker$cluster)
            top_1_markers <- lapply(top_1_markers, function(x){
                x <- x[order(x[,wt], decreasing = T),]
                x <- x[1,]
            })
            top_1_markers <- do.call(rbind.data.frame, top_1_markers)
            
            Idents(data) <- "seurat_clusters"
            results[['p61plots']] <- RidgePlot(data, features = unique(unlist(top_1_markers$gene)), ncol = 4,
                                               cols = cluster_colors)+
                plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
            
            current <- data@meta.data
            total_counts <- data.frame(table(current[,integration_cluster]))
            
            current <- data.frame(table(current[,c("SID",integration_cluster)]))
            current <- current[current$Freq != 0,]
            colnames(current) <- c("SID","CLUSTER","COUNT")
            current$CLUSTER_COUNT <- total_counts[match(current$CLUSTER, total_counts$Var1),"Freq"]
            current$PROPORTION <- current$COUNT/current$CLUSTER_COUNT
            
            node_proportion <- current
            
            results[['p62plots']] <- ggplot(node_proportion, aes(CLUSTER, PROPORTION, fill = SID))+
                geom_bar(stat="identity", alpha=0.8)+
                coord_polar()+
                scale_fill_viridis(discrete = T, option = "C")+
                ggtitle(paste("Frequency of Samples in Each Node\n", project_name, sep = ""))+
                theme_bw(base_size = 28)+
                theme(plot.margin = unit(c(3,3,3,3), "cm"),
                      plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                      legend.title=element_text(size=20, face = "bold"),
                      legend.text=element_text(size=15),
                      legend.key.size = unit(2, 'lines'),
                      axis.title.x = element_text(colour="black", size = 20, face = "bold", vjust = -10),
                      axis.title.y = element_text(colour="black", size = 20, face = "bold", vjust = 10),
                      strip.text = element_text(size = 20, face = "bold"),
                      axis.text.x=element_text(colour="black", size = 20),
                      axis.text.y=element_text(colour="black", size = 20),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
            
            current <- data@meta.data
            total_counts <- data.frame(table(current$CELL_TYPE))
            
            current <- data.frame(table(current[,c("SID","CELL_TYPE")]))
            current <- current[current$Freq != 0,]
            colnames(current) <- c("SID","CELL_TYPE","COUNT")
            current$CELL_TYPE_COUNT <- total_counts[match(current$CELL_TYPE, total_counts$Var1),"Freq"]
            current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT

            results[['p63plots']] <- ggplot(current,aes(SID,PROPORTION,fill=CELL_TYPE))+
                geom_bar(stat = "identity", position = "fill", color = "black", size = 1)+ #, color = "black", size = 1
                theme_classic()+
                # facet_wrap(~GROUP)+
                theme(plot.margin = unit(c(1,1,1,1), "cm"),
                      axis.text.x = element_text(angle = 45, size = 15, hjust = 1, vjust = 1),
                      axis.text.y = element_text(size = 20),
                      axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
                      axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
                      legend.title = element_text(size = 20),
                      legend.text = element_text(size = 15),
                      legend.key.size = unit(1, "cm"),
                      legend.position = "right",
                      strip.text.x = element_text(size = 25),
                      strip.text.y = element_text(size = 25),
                      plot.title = element_text(size = 20, face=1, hjust = 0.5))+
                scale_fill_manual(values = ct_colors)+
                guides(color=guide_legend(title="CELL TYPES", ncol = 1),fill = guide_legend(title="CELL TYPES", ncol = 1))
            
            #####################################################################################################################
            current <- group_medianexpr(current_out$current_data_markers, data, group = "CELL_TYPE", cell_type = T)
            plot_median_cell_type <- current$plot_median
            top_markers_cell_type <- current$top_markers
            plot_median_cell_type[which(is.na(plot_median_cell_type))] <- 0
            results[['p64plots']] <- complex_heatmap(plot_median_cell_type, col = color_conditions$BlueYellowRed)
            
            p1 <- NULL
            p2 <- NULL
            p3 <- NULL
            cell_types <- sort(as.character(unique(top_markers_cell_type$CELL_TYPE)))
            
            plotx <- gen10x_plotx(data, groups = NULL)
            plotx$CELL_TYPE <- data$CELL_TYPE
            
            cct_colors <- gen_colors(color_conditions$tenx,length(unique(plotx$CELL_TYPE)))
            names(cct_colors) <- sort(unique(plotx$CELL_TYPE))
            
            p1 <- ggplot(data=data.frame(top_markers_cell_type), aes_string(x="gene", y=wt, fill="CELL_TYPE")) +
                geom_bar(stat="identity", position=position_dodge())+
                scale_fill_manual(values = cct_colors)+
                theme_classic() +
                coord_flip()+
                ylab(wt)+
                facet_wrap(~CELL_TYPE, scales = "free_y") +
                theme(axis.text.y=element_text(size = 7))
            
            p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CELL_TYPE",
                               plot_title = integration_name,col = ct_colors,numeric = F,
                               annot = T, legend_position = "right", point_size = 1)
            p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CELL_TYPE",
                               plot_title = integration_name,col = ct_colors,numeric = F,
                               annot = T, legend_position = "right", point_size = 1)
            
            results[['p65plots']] <- p1
            results[['p652plots']] <- p2
            results[['p653plots']] <- p3
            
            # grid.arrange(p1, p2, ncol = 2)
            # grid.arrange(p1, p3, ncol = 2)
            
            current_out$current_data_markers$CELL_TYPE <- data@meta.data[match(current_out$current_data_markers$cluster, data@meta.data[,integration_cluster]),"CELL_TYPE"]
            # top_1_markers <- current_out$current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = 1, eval(parse(text = wt)))
            
            top_1_markers <- split(current_out$current_data_marker, current_out$current_data_marker$CELL_TYPE)
            top_1_markers <- lapply(top_1_markers, function(x){
                x <- x[order(x[,wt], decreasing = T),]
                x <- x[1,]
            })
            top_1_markers <- do.call(rbind.data.frame, top_1_markers)
            
            Idents(data) <- "CELL_TYPE"
            results[['p66plots']] <- RidgePlot(data, features = unique(unlist(top_1_markers$gene)), ncol = 2,
                                               cols = ct_colors)+
                plot_annotation(title = integration_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
            
            current <- data@meta.data
            total_counts <- data.frame(table(current$CELL_TYPE))
            
            current <- data.frame(table(current[,c("SID","CELL_TYPE")]))
            current <- current[current$Freq != 0,]
            colnames(current) <- c("SID","CELL_TYPE","COUNT")
            current$CELL_TYPE_COUNT <- total_counts[match(current$CELL_TYPE, total_counts$Var1),"Freq"]
            current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT
            
            node_proportion <- current
            
            results[['p67plots']] <- ggplot(node_proportion, aes(CELL_TYPE, PROPORTION, fill = SID))+
                geom_bar(stat="identity", alpha=0.8)+
                coord_polar()+
                scale_fill_viridis(discrete = T)+
                ggtitle(paste("Frequency of Samples in Each Cell Type\n", project_name, sep = ""))+
                theme_bw(base_size = 25)+
                theme(plot.margin = unit(c(3,3,3,3), "cm"),
                      plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                      legend.title=element_text(size=20, face = "bold"),
                      legend.text=element_text(size=15),
                      legend.key.size = unit(2, 'lines'),
                      axis.title.x = element_text(colour="black", size = 20, face = "bold", vjust = -10),
                      axis.title.y = element_text(colour="black", size = 20, face = "bold", vjust = 10),
                      strip.text = element_text(size = 20, face = "bold"),
                      axis.text.x=element_text(colour="black", size = 20),
                      axis.text.y=element_text(colour="black", size = 20),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
            
            Idents(data) <- integration_cluster
            
            ############ GROUP COMPARISON #########################################################
            
            Idents(data) <- integration_cluster
            current_clusters <- sort(as.numeric(as.character(unique(Idents(data)))), decreasing = F)
            
            if(length(unique(data$GROUP)) > 1){
                current_out <- deanalysis(data, current_clusters, plot_title = project_name,group="GROUP",de_analysis = "finddegroup")
                results[['p68data']] <- current_out$current_data_markers
                de_type <- current_out$de_type
                de_name <- current_out$de_name
                top1 <- current_out$top1
                topn <- current_out$topn
                current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
                
                results[['p69plots']] <- DotPlot(data, features = unique(top1$gene), cols = cgroup_colors, dot.scale = 8,
                                                 split.by = "GROUP", cluster.idents = F) + RotatedAxis()+ylab("CLUSTER+GROUP")+
                    plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5)))
                
                for(k in 1:length(current_clusters)){
                    current <- topn[which(topn$cluster == current_clusters[k]),]
                    current <- current[order(current$p_val_adj, decreasing = F),]
                    results[['p70plots']][[length(results[['p70plots']])+1]] <- VlnPlot(data, features = unique(current$gene), pt.size = 0, split.by = "GROUP",stack = T,flip = T,
                                                                                        cols = cgroup_colors)&
                        xlab("CLUSTERS")&
                        plot_annotation(title = paste(project_name, ": CLUSTER ", current_clusters[k], sep = ""),
                                        theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
                    names(results[['p70plots']])[length(results[['p70plots']])] <-paste("CLUSTER_",current_clusters[k], sep = "")
                }
                
                for(i in 1:length(current_out$featureplot)){
                    results[['p71plots']][[length(results[['p71plots']])+1]] <- current_out$featureplot[[i]]
                    names(results[['p71plots']])[length(results[['p71plots']])] <- names(current_out$featureplot)[i]
                    
                }
                
                # current_clusters <- sort(as.numeric(as.character(unique(Idents(data)))), decreasing = F)
                # Idents(data) <- integration_cluster
                # current_out <- deanalysis(data, current_clusters, plot_title = project_name,group="GROUP",de_analysis = "findconservedmarkers")
                # results[['p72data']] <- current_out$current_data_markers
                # de_type <- current_out$de_type
                # de_name <- current_out$de_name
                # top1 <- current_out$top1
                # topn <- current_out$topn
                # current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
                #
                # for(i in 1:length(current_out$featureplot)){
                #   results[['p73plots']][[length(results[['p73plots']])+1]] <- current_out$featureplot[[i]]
                #   names(results[['p73plots']])[length(results[['p73plots']])] <- paste("CLUSTER_",names(current_out$featureplot)[i], sep = "")
                # }
                
            }
        }
        
        #################### PSEUDOTIME TRAJECTORY #######################################################
        mono3_current <- NULL
        sce_current <- NULL
        
        for(j in 1:length(data_current)){
            mono3_current[[j]] <- as.cell_data_set(data_current[[j]])
            
            ##### UMAP ######
            reducedDim(mono3_current[[j]], type = "PCA") <- data_current[[j]]@reductions$pca@cell.embeddings
            mono3_current[[j]]@reduce_dim_aux$prop_var_expl <- data_current[[j]]@reductions$pca@stdev
            mono3_current[[j]]@reduce_dim_aux$gene_loadings <- data_current[[j]]@reductions[["pca"]]@feature.loadings
            mono3_current[[j]]@int_colData@listData$reducedDims$UMAP <- data_current[[j]]@reductions$umap@cell.embeddings
            mono3_current[[j]]@clusters$UMAP$clusters <- data_current[[j]]$seurat_clusters
            mono3_current[[j]]@clusters@listData[["UMAP"]][["clusters"]] <- data_current[[j]]$seurat_clusters
            mono3_current[[j]]@clusters@listData[["UMAP"]]$cluster_result$optim_res$membership <- data_current[[j]]$seurat_clusters
            mono3_current[[j]]@clusters@listData[["UMAP"]][["louvain_res"]] <- NULL
            rownames(mono3_current[[j]]@principal_graph_aux$UMAP$dp_mst) <- NULL
            colnames(mono3_current[[j]]@int_colData@listData$reducedDims$UMAP) <- NULL
            mono3_current[[j]] <- cluster_cells(mono3_current[[j]], reduction_method = "UMAP", resolution = 1e-3)
            mono3_current[[j]]@clusters@listData[["UMAP"]][["clusters"]] <- data_current[[j]]$seurat_clusters
            mono3_current[[j]]@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
            mono3_current[[j]] <- learn_graph(mono3_current[[j]], use_partition = F)
            
            ccolors <- gen_colors(color_conditions$colorful,length(unique(mono3_current[[j]]@clusters$UMAP$clusters)))
            names(ccolors) <- sort(as.numeric(as.character(unique(mono3_current[[j]]@clusters$UMAP$clusters))))
            
            results[['p72plots']][[j]] <- plot_pseudo(mono3_current[[j]], reduction = "UMAP",
                                                      group = "seurat_clusters", label_size = 7,
                                                      paste(annot_names[j],"\nTRAJECTORY - COLOR BY CLUSTERS"),
                                                      ccolors, length(ccolors))
            names(results[['p72plots']])[j] <- names(data_current)[j]
            
            results[['p73plots']][[j]] <- ggplotly(results[['p72plots']][[j]]+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                                       guides(color=guide_legend(title="CLUSTER"))+
                                                       ggtitle(paste("<b>",names(data_current)[j],"</b>",sep = "")))
            names(results[['p73plots']])[j] <- names(data_current)[j]
            
            results[['p74plots']][[j]] <- plot_pseudo(mono3_current[[j]], reduction = "UMAP",
                                                      group = "CELL_TYPE", label_size = 6,
                                                      paste(annot_names[j],"\nTRAJECTORY - COLOR BY CELL TYPE"),
                                                      color_conditions$tenx, length(unique(mono3_current[[j]]$CELL_TYPE)))
            names(results[['p74plots']])[j] <- names(data_current)[j]
            
            results[['p75plots']][[j]] <- ggplotly(results[['p74plots']][[j]]+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                                       # guides(color=guide_legend(title=""))+
                                                       ggtitle(paste("<b>",names(data_current)[j],"</b>",sep = "")))
            names(results[['p75plots']])[j] <- names(data_current)[j]
            
            results[['p76plots']][[j]] <- plot_cells(mono3_current[[j]],
                                                     color_cells_by = "partition",
                                                     group_label_size = 7,
                                                     graph_label_size = 3,
                                                     cell_size = 1,
                                                     # cell_stroke = I(2/2),
                                                     alpha = 1,
                                                     trajectory_graph_segment_size = 1.5,
                                                     label_cell_groups=FALSE,
                                                     label_leaves=TRUE,
                                                     label_branch_points=TRUE)
            
            results[['p76plots']][[j]] <- adjust_plot(results[['p76plots']][[j]],col = color_conditions$bright, n = length(unique(mono3_current[[j]]@clusters$UMAP$clusters)),
                                                      plot_title = paste(annot_names[j],"\nTRAJECTORY - COLOR BY PARTITION"))
            
            names(results[['p76plots']])[j] <- names(data_current)[j]
            
            current_clusters <- sort(as.numeric(as.character(unique(mono3_current[[j]]@clusters$UMAP$clusters))), decreasing = F)
            for(i in 1:length(current_clusters)){
                mono3_current[[j]] <- order_cells(mono3_current[[j]],
                                                  root_pr_nodes=get_earliest_principal_node(mono3_current[[j]],
                                                                                            group = "seurat_clusters", group_element = current_clusters[i]))
                
                results[['p77plots']][[length(results[['p77plots']])+1]] <- plot_cells(mono3_current[[j]],
                                                                                       color_cells_by = "pseudotime",
                                                                                       group_label_size = 7,
                                                                                       graph_label_size = 5,
                                                                                       cell_size = 1,
                                                                                       # cell_stroke = I(2/2),
                                                                                       alpha = 1,
                                                                                       trajectory_graph_segment_size = 1.5,
                                                                                       label_cell_groups=FALSE,
                                                                                       label_leaves=FALSE,
                                                                                       label_branch_points=FALSE)
                
                results[['p77plots']][[length(results[['p77plots']])]] <- adjust_plot(results[['p77plots']][[length(results[['p77plots']])]], col = color_conditions$colorful, n = length(current_clusters),
                                                                                      plot_title = paste(annot_names[j],"\nPSEUDOTIME - ROOT: CLUSTER ", current_clusters[i], sep = ""),
                                                                                      fill = T)
                
                names(results[['p77plots']])[length(results[['p77plots']])] <- paste("CLUSTER_",current_clusters[i],"|",names(data_current)[j], sep = "")
                
            }
            
            pseudo_para <- c("PSEUDOTIME_PC1","PSEUDOTIME_UMAP1","PSEUDOTIME_tSNE1")
            
            # if(to_annotate[[j]] == TRUE){
            chosen_para <- c("seurat_clusters", "CELL_TYPE")
            chosen_para_names <- c("CLUSTERS", "CELL_TYPE")
            
            plotx <- data.frame(
                PSEUDOTIME_PC1 = rank(data_current[[j]]@reductions$pca@cell.embeddings[,"PC_1"]),
                PSEUDOTIME_UMAP1 = rank(data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_1"]),
                PSEUDOTIME_tSNE1 = rank(data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_1"]),
                data_current[[j]]@meta.data)
            
            k <- 2
            for(m in 1:length(chosen_para)){
                for(n in 1:length(pseudo_para)){
                  cp <- NULL
                  pseudocols <- NULL
                  pseudocols <- prepare_cols(if(chosen_para_names[m] == "CELL_TYPE"){color_conditions$tenx}else{color_conditions$colorful}, length(unique(plotx[,chosen_para[m]])))
                  plotx[,chosen_para[m]] <- factor(plotx[,chosen_para[m]], levels = unique(unlist(plotx[order(unlist(plotx[,pseudo_para[n]]), decreasing = F),chosen_para[m]])))
                  cp <- ggplot(plotx, aes_string(x = pseudo_para[n], y = chosen_para[m], colour = chosen_para[m])) + geom_quasirandom(groupOnX = FALSE) +
                    scale_color_manual(values = pseudocols) +
                    ylab(chosen_para_names[m])+
                    ggtitle(annot_names[j])
                  cp <- adjust_theme(cp, legend = "none")
                  results[['p78plots']][[length(results[['p78plots']])+1]] <- cp
                    
                    names(results[['p78plots']])[length(results[['p78plots']])] <- paste(names(data_current)[j],"|",chosen_para_names[m],"|",pseudo_para[n], sep = "")
                    
                    k <- k + 1
                }
            }
            
            # if(to_annotate[[j]] == TRUE){
            current_cell_types <- sort(as.character(unique(colData(mono3_current[[j]])[,"CELL_TYPE"])))
            for(i in 1:length(current_cell_types)){
                mono3_current[[j]] <- order_cells(mono3_current[[j]],
                                                  root_pr_nodes=get_earliest_principal_node(mono3_current[[j]],
                                                                                            group = "CELL_TYPE", group_element = current_cell_types[i]))
                
                results[['p79plots']][[length(results[['p79plots']])+1]] <- plot_cells(mono3_current[[j]],
                                                                                       color_cells_by = "pseudotime",
                                                                                       group_label_size = 7,
                                                                                       graph_label_size = 5,
                                                                                       cell_size = 1,
                                                                                       # cell_stroke = I(2/2),
                                                                                       alpha = 1,
                                                                                       trajectory_graph_segment_size = 1.5,
                                                                                       label_cell_groups=FALSE,
                                                                                       label_leaves=FALSE,
                                                                                       label_branch_points=FALSE)
                
                results[['p79plots']][[length(results[['p79plots']])]] <- adjust_plot(results[['p79plots']][[length(results[['p79plots']])]], col = color_conditions$tenx, n = length(current_cell_types),
                                                                                      plot_title = paste(annot_names[j],"\nPSEUDOTIME - ROOT CELL TYPE: ", toupper(current_cell_types[i]), sep = ""),
                                                                                      fill = T)
                names(results[['p79plots']])[length(results[['p79plots']])] <- paste(names(data_current)[j],"|",toupper(current_cell_types[i]), sep = "")
            }
            
            data_current[[j]] <- RunUMAP(data_current[[j]], n.components = 3, reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
            
            results[['p80plots']][[j]] <- plot_3d(data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                                                  data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_2"],
                                                  data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_3"],
                                                  color_group = data_current[[j]]@meta.data[,"CELL_TYPE"],
                                                  plot_title = names(data_current)[j],
                                                  col = color_conditions$tenx,
                                                  n = length(unique(data_current[[j]]$seurat_clusters)),
                                                  lt = "CELL TYPE", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                                  current_text= paste("Cluster ",data_current[[j]]$seurat_clusters,"\n",
                                                                      "Cell: ",row.names(data_current[[j]]@meta.data),"\n",
                                                                      data_current[[j]]@meta.data[,"CELL_TYPE"], sep = ""))
            names(results[['p80plots']])[j] <- names(data_current)[j]
            
            data_current[[j]] <- RunTSNE(data_current[[j]], dim.embed = 3, reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30), check_duplicates = FALSE)
            
            results[['p81plots']][[j]] <- plot_3d(data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                                  data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                                  data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                                  color_group = data_current[[j]]@meta.data[,"CELL_TYPE"],
                                                  plot_title = names(data_current)[j],
                                                  col = color_conditions$tenx,
                                                  n = length(unique(data_current[[j]]$seurat_clusters)),
                                                  lt = "CELL TYPE", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                                  current_text= paste("Cluster ",data_current[[j]]$seurat_clusters,"\n",
                                                                      "Cell: ",row.names(data_current[[j]]@meta.data),"\n",
                                                                      data_current[[j]]@meta.data[,"CELL_TYPE"], sep = ""))
            names(results[['p81plots']])[j] <- names(data_current)[j]
            
            results[['p82plots']][[j]] <- plot_3d(data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                                                  data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_2"],
                                                  data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_3"],
                                                  color_group = data_current[[j]]$seurat_clusters,
                                                  plot_title = names(data_current)[j],
                                                  col = color_conditions$colorful,
                                                  n = length(unique(data_current[[j]]$seurat_clusters)),
                                                  lt = "CLUSTER", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                                  current_text= paste("Cell: ",row.names(data_current[[j]]@meta.data),"\nCluster: ",
                                                                      data_current[[j]]$seurat_clusters, sep = ""))
            names(results[['p82plots']])[j] <- names(data_current)[j]
            
            results[['p83plots']][[j]] <- plot_3d(data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                                  data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                                  data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                                  color_group = data_current[[j]]$seurat_clusters,
                                                  plot_title = names(data_current)[j],
                                                  col = color_conditions$colorful,
                                                  n = length(unique(data_current[[j]]$seurat_clusters)),
                                                  lt = "CLUSTER", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                                  current_text= paste("Cell: ",row.names(data_current[[j]]@meta.data),"\nCluster: ",
                                                                      data_current[[j]]$seurat_clusters, sep = ""))
            names(results[['p83plots']])[j] <- names(data_current)[j]
            
            data_current[[j]] <- RunTSNE(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30), check_duplicates = FALSE)
            data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:ifelse(length(data_current[[j]]@reductions$pca) < 30, length(data_current[[j]]@reductions$pca), 30))
        }
        
        if(length(data_current) > 1){
            
            pseudo_para <- c("PSEUDOTIME_PC1","PSEUDOTIME_UMAP1","PSEUDOTIME_tSNE1")
            chosen_para <- c(integration_cluster, "CELL_TYPE","orig.ident")
            chosen_para_names <- c("CLUSTERS", "CELL_TYPE","SAMPLE_ID")
            plotx <- data.frame(
                PSEUDOTIME_PC1 = rank(data@reductions$pca@cell.embeddings[,"PC_1"]),
                PSEUDOTIME_UMAP1 = rank(data@reductions$umap@cell.embeddings[,"UMAP_1"]),
                PSEUDOTIME_tSNE1 = rank(data@reductions$tsne@cell.embeddings[,"tSNE_1"]),
                data@meta.data)
            
            for(m in 1:length(chosen_para)){
                for(n in 1:length(pseudo_para)){
                  cp <- NULL
                  pseudocols <- NULL
                  pseudocols <- prepare_cols(col = NULL, length(unique(plotx[,chosen_para[m]])))
                  cp <- ggplot(plotx, aes_string(x = pseudo_para[n], y = chosen_para[m], colour = chosen_para[m])) + geom_quasirandom(groupOnX = FALSE) +
                    scale_color_manual(values = pseudocols) +
                    ylab(chosen_para_names[m])+
                    ggtitle(project_name)
                  cp <- adjust_theme(cp, legend = "none")
                  
                    results[['p84plots']][[length(results[['p84plots']])+1]] <- cp
                    names(results[['p84plots']])[length(results[['p84plots']])] <- paste(pseudo_para[n], "|", chosen_para_names[m],sep = "")
                }
            }
            
            if(integration_name == "SEURAT_INTEGRATED"){
                DefaultAssay(data) <- "integrated"
                mono3data <- as.cell_data_set(data)
            }else if(integration_name == "HARMONY_INTEGRATED"){
                DefaultAssay(data) <- "RNA"
                mono3data <- as.cell_data_set(data)
            }
            
            ##### UMAP ######
            named_clusters <- data@meta.data[,integration_cluster]
            names(named_clusters) <- row.names(data@meta.data)
            reducedDim(mono3data, type = "PCA") <- data@reductions[[reduction_method]]@cell.embeddings
            mono3data@reduce_dim_aux$prop_var_expl <- data@reductions$pca@stdev
            mono3data@reduce_dim_aux$gene_loadings <- data@reductions[[reduction_method]]@feature.loadings
            mono3data@int_colData@listData$reducedDims$UMAP <- data@reductions$umap@cell.embeddings
            mono3data@clusters$UMAP$clusters <- named_clusters
            mono3data@clusters@listData[["UMAP"]][["clusters"]] <- named_clusters
            mono3data@clusters@listData[["UMAP"]]$cluster_result$optim_res$membership <- named_clusters
            mono3data@clusters@listData[["UMAP"]][["louvain_res"]] <- NULL
            rownames(mono3data@principal_graph_aux$UMAP$dp_mst) <- NULL
            colnames(mono3data@int_colData@listData$reducedDims$UMAP) <- NULL
            mono3data <- cluster_cells(mono3data, reduction_method = "UMAP", resolution = 1e-3)
            mono3data@clusters@listData[["UMAP"]][["clusters"]] <- named_clusters
            mono3data@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
            mono3data <- learn_graph(mono3data, use_partition = F)
            
            results[['p85plots']] <- plot_pseudo(mono3data, reduction = "UMAP",
                                                 group = integration_cluster, label_size = 7,
                                                 # cell_size = 0.5,traj_size = 0.8,
                                                 paste(project_name,"\nINTEGRATED TRAJECTORY - COLOR BY CLUSTERS"),
                                                 color_conditions$colorful, length(unique(mono3data@clusters$UMAP$clusters)))+facet_wrap(~GROUP)
            
            results[['p86plots']] <- ggplotly(results[['p85plots']]+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                                  guides(color=guide_legend(title="CLUSTER"))+
                                                  ggtitle(paste("<b>",project_name,"</b>",sep = "")))
            
            results[['p87plots']] <- plot_pseudo(mono3data, reduction = "UMAP",
                                                 group = "CELL_TYPE", label_size = 6,
                                                 # cell_size = 0.5,traj_size = 0.8,
                                                 paste(project_name,"\nINTEGRATED TRAJECTORY - COLOR BY CELL TYPE"),
                                                 color_conditions$tenx, length(unique(mono3data$CELL_TYPE)))+
                facet_wrap(~GROUP)
            
            results[['p88plots']] <- ggplotly(results[['p87plots']]+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                                  # guides(color=guide_legend(title=""))+
                                                  ggtitle(paste("<b>",project_name,"</b>",sep = "")))
            
            results[['p89plots']] <- plot_cells(mono3data,
                                                color_cells_by = "partition",
                                                group_label_size = 7,
                                                graph_label_size = 3,
                                                cell_size = 1,
                                                # cell_stroke = I(2/2),
                                                alpha = 1,
                                                trajectory_graph_segment_size = 1.5,
                                                label_cell_groups=FALSE,
                                                label_leaves=TRUE,
                                                label_branch_points=TRUE)
            
            results[['p89plots']] <- adjust_plot(results[['p89plots']],col = color_conditions$colorful, n = length(unique(mono3data@clusters$UMAP$clusters)),
                                                 plot_title = paste(project_name,"\nTRAJECTORY - COLOR BY PARTITION"))
            
            current_clusters <- sort(as.numeric(as.character(unique(mono3data@clusters$UMAP$clusters))), decreasing = F)
            for(i in 1:length(current_clusters)){
                mono3data <- order_cells(mono3data,
                                         root_pr_nodes=get_earliest_principal_node(mono3data,
                                                                                   group = "seurat_clusters", group_element = current_clusters[i]))
                
                results[['p90plots']][[length(results[['p90plots']])+1]] <- plot_cells(mono3data,
                                                                                       color_cells_by = "pseudotime",
                                                                                       group_label_size = 7,
                                                                                       graph_label_size = 5,
                                                                                       cell_size = 1,
                                                                                       # cell_stroke = I(2/2),
                                                                                       alpha = 1,
                                                                                       trajectory_graph_segment_size = 1.5,
                                                                                       label_cell_groups=FALSE,
                                                                                       label_leaves=FALSE,
                                                                                       label_branch_points=FALSE)+
                    facet_wrap(~GROUP)
                
                results[['p90plots']][[length(results[['p90plots']])]] <- adjust_plot(results[['p90plots']][[length(results[['p90plots']])]], col = color_conditions$colorful, n = length(current_clusters),
                                                                                      plot_title = paste(project_name,"\nPSEUDOTIME - ROOT: CLUSTER ", current_clusters[i], sep = ""),
                                                                                      fill = T)
                names(results[['p90plots']])[length(results[['p90plots']])] <- current_clusters[i]
            }
            
            if(reduction_method == "harmony"){
                DefaultAssay(data) <- 'RNA'
                data <- RunTSNE(data, dim.embed = 3, reduction = "harmony", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30), check_duplicates = FALSE)
                data <- RunUMAP(data, n.components = 3, reduction = "harmony", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30))
                
            }else if(reduction_method == "pca"){
                DefaultAssay(data) <- 'integrated'
                data <- RunTSNE(data, dim.embed = 3, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30), check_duplicates = FALSE)
                data <- RunUMAP(data, n.components = 3, reduction = "pca", dims = 1:ifelse(length(data@reductions$pca) < 30, length(data@reductions$pca), 30))
            }
            
            results[['p91plots']] <- plot_3d(data@reductions$umap@cell.embeddings[,"UMAP_1"],
                                             data@reductions$umap@cell.embeddings[,"UMAP_2"],
                                             data@reductions$umap@cell.embeddings[,"UMAP_3"],
                                             color_group = data@meta.data[,"CELL_TYPE"],
                                             plot_title = names(data_current)[j],
                                             col = color_conditions$tenx,
                                             n = length(unique(data$CELL_TYPE)),
                                             lt = "CELL TYPE", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                             current_text= paste("Cluster ",data@meta.data[,integration_cluster],"\n",
                                                                 "Cell: ",row.names(data@meta.data),"\n",
                                                                 data@meta.data[,"CELL_TYPE"], sep = ""))
            
            results[['p92plots']] <- plot_3d(data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                             data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                             data@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                             color_group = data@meta.data[,"CELL_TYPE"],
                                             plot_title = names(data_current)[j],
                                             col = color_conditions$tenx,
                                             n = length(unique(data$CELL_TYPE)),
                                             lt = "CELL TYPE", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                             current_text= paste("Cluster ",data@meta.data[,integration_cluster],"\n",
                                                                 "Cell: ",row.names(data@meta.data),"\n",
                                                                 data@meta.data[,"CELL_TYPE"], sep = ""))
            
            
            results[['p93plots']] <- plot_3d(data@reductions$umap@cell.embeddings[,"UMAP_1"],
                                             data@reductions$umap@cell.embeddings[,"UMAP_2"],
                                             data@reductions$umap@cell.embeddings[,"UMAP_3"],
                                             color_group = data@meta.data[,integration_cluster],
                                             plot_title = names(data_current)[j],
                                             col = color_conditions$colorful,
                                             n = length(unique(data@meta.data[,integration_cluster])),
                                             lt = "CLUSTER", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                             current_text= paste("Cell: ",row.names(data@meta.data),"\nCluster: ",
                                                                 data@meta.data[,integration_cluster], sep = ""))
            
            results[['p94plots']] <- plot_3d(data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                             data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                             data@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                             color_group = data@meta.data[,integration_cluster],
                                             plot_title = names(data_current)[j],
                                             col = color_conditions$colorful,
                                             n = length(unique(data@meta.data[,integration_cluster])),
                                             lt = "CLUSTER", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                             current_text= paste("Cell: ",row.names(data@meta.data),"\nCluster: ",
                                                                 data@meta.data[,integration_cluster], sep = ""))
        }
        ##########################################################################################
        #######################################################################################################################################
        # saveRDS(results,"results.RDS")
        updateSelectInput(session, inputId = 'p1id', label = 'Choose a sample to display', choices = names(results[['p1plots']]), selected = names(results[['p1plots']])[1])
        updateSelectInput(session, inputId = 'p2id', label = 'Choose a sample to display', choices = names(results[['p2plots']]), selected = names(results[['p2plots']])[1])
        updateSelectInput(session, inputId = 'p3id', label = 'Choose a sample to display', choices = names(results[['p3plots']][which(!is.na(names(results[['p3plots']])))]), selected = names(results[['p3plots']][which(!is.na(names(results[['p3plots']])))])[1])
        updateSelectInput(session, inputId = 'p4id', label = 'Choose a sample to display', choices = names(results[['p4plots']][which(!is.na(names(results[['p4plots']])))]), selected = names(results[['p4plots']][which(!is.na(names(results[['p4plots']])))])[1])
        updateSelectInput(session, inputId = 'p5id', label = 'Choose a sample to display', choices = names(results[['p5plots']][which(!is.na(names(results[['p5plots']])))]), selected = names(results[['p5plots']][which(!is.na(names(results[['p5plots']])))])[1])
        updateSelectInput(session, inputId = 'p6id', label = 'Choose a sample to display', choices = names(results[['p6plots']][which(!is.na(names(results[['p6plots']])))]), selected = names(results[['p6plots']][which(!is.na(names(results[['p6plots']])))])[1])
        updateSelectInput(session, inputId = 'p7id', label = 'Choose a sample to display', choices = names(results[['p7plots']][which(!is.na(names(results[['p7plots']])))]), selected = names(results[['p7plots']][which(!is.na(names(results[['p7plots']])))])[1])
        updateSelectInput(session, inputId = 'p8id', label = 'Choose a sample to display', choices = names(results[['p8plots']]), selected = names(results[['p8plots']])[1])
        updateSelectInput(session, inputId = 'p9id', label = 'Choose a sample to display', choices = names(results[['p9plots']]), selected = names(results[['p9plots']])[1])
        # updateSelectInput(session, inputId = 'p10id', label = 'Choose a sample to display', choices = names(results[['p10plots']]), selected = names(results[['p10plots']])[1])
        updateSelectInput(session, inputId = 'p11id', label = 'Choose a sample to display', choices = names(results[['p11plots']]), selected = names(results[['p11plots']])[1])
        updateSelectInput(session, inputId = 'p12id', label = 'Choose a sample to display', choices = names(results[['p12plots']]), selected = names(results[['p12plots']])[1])
        updateSelectInput(session, inputId = 'p13id', label = 'Choose a sample to display', choices = names(results[['p13plots']]), selected = names(results[['p13plots']])[1])
        updateSelectInput(session, inputId = 'p14id', label = 'Choose a sample to display', choices = names(results[['p14plots']]), selected = names(results[['p14plots']])[1])
        updateSelectInput(session, inputId = 'p15id', label = 'Choose a sample to display', choices = names(results[['p15plots']]), selected = names(results[['p15plots']])[1])
        updateSelectInput(session, inputId = 'p16id', label = 'Choose a sample to display', choices = names(results[['p16plots']]), selected = names(results[['p16plots']])[1])
        updateSelectInput(session, inputId = 'p17id', label = 'Choose a sample to display', choices = names(results[['p17data']]), selected = names(results[['p17data']])[1])
        updateSelectInput(session, inputId = 'p18id', label = 'Choose a sample to display', choices = names(results[['p18plots']]), selected = names(results[['p18plots']])[1])
        updateSelectInput(session, inputId = 'p19id', label = 'Choose a sample to display', choices = names(results[['p19plots']]), selected = names(results[['p19plots']])[1])
        updateSelectInput(session, inputId = 'p193id', label = 'Choose a sample to display', choices = names(results[['p193plots']]), selected = names(results[['p193plots']])[1])
        updateSelectInput(session, inputId = 'p20id', label = 'Choose a sample to display', choices = names(results[['p20plots']]), selected = names(results[['p20plots']])[1])
        updateSelectInput(session, inputId = 'p21id', label = 'Choose a sample to display', choices = names(results[['p21plots']]), selected = names(results[['p21plots']])[1])
        p22samples <- gsub(".*\\|(.*)","\\1",names(results[['p22plots']]))
        p22clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p22plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p22id', label = 'Choose a sample to display', choices = p22samples, selected = p22samples[1])
        updateSliderInput(session, inputId ="p222id", label = 'Choose a cluster to display', value = 0, min = min(p22clusters), max = max(p22clusters), step = 1)
        updateSelectInput(session, inputId = 'p23id', label = 'Choose a sample to display', choices = names(results[['p23plots']]), selected = names(results[['p23plots']])[1])
        updateSelectInput(session, inputId = 'p24id', label = 'Choose a sample to display', choices = names(results[['p24plots']]), selected = names(results[['p24plots']])[1])
        updateSelectInput(session, inputId = 'p25id', label = 'Choose a sample to display', choices = names(results[['p25plots']]), selected = names(results[['p25plots']])[1])
        p26samples <- gsub(".*\\|(.*)","\\1",names(results[['p26plots']]))
        p262clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p26plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p26id', label = 'Choose a sample to display', choices = p26samples, selected = p26samples[1])
        updateSliderInput(session, inputId ="p262id", label = 'Choose a cluster to display', value = 0, min = min(p262clusters), max = max(p262clusters), step = 1)
        updateSelectInput(session, inputId = 'p27id', label = 'Choose a sample to display', choices = names(results[['p27data']]), selected = names(results[['p27data']])[1])
        p28samples <- gsub(".*\\|(.*)","\\1",names(results[['p28plots']]))
        p282clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p28plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p28id', label = 'Choose a sample to display', choices = p28samples, selected = p28samples[1])
        updateSliderInput(session, inputId ="p282id", label = 'Choose a cluster to display', value = 0, min = min(p282clusters), max = max(p282clusters), step = 1)
        p29samples <- gsub(".*\\|(.*)","\\1",names(results[['p29plots']]))
        p292clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p29plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p29id', label = 'Choose a sample to display', choices = p29samples, selected = p29samples[1])
        updateSliderInput(session, inputId ="p292id", label = 'Choose a cluster to display', value = 0, min = min(p292clusters), max = max(p292clusters), step = 1)
        p30samples <- gsub(".*\\|(.*)","\\1",names(results[['p30plots']]))
        p302clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p30plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p30id', label = 'Choose a sample to display', choices = p30samples, selected = p30samples[1])
        updateSliderInput(session, inputId ="p302id", label = 'Choose a cluster to display', value = 0, min = min(p302clusters), max = max(p302clusters), step = 1)
        p31samples <- unique(gsub(".*\\|(.*)","\\1",names(results[['p31plots']])))
        p312clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p31plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p31id', label = 'Choose a sample to display', choices = p31samples, selected = p31samples[1])
        updateSliderInput(session, inputId ="p312id", label = 'Choose a cluster to display', value = 0, min = min(p312clusters), max = max(p312clusters), step = 1)
        p32samples <- unique(gsub(".*\\|(.*)","\\1",names(results[['p32plots']])))
        p322clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p32plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p32id', label = 'Choose a sample to display', choices = p32samples, selected = p32samples[1])
        updateSliderInput(session, inputId ="p322id", label = 'Choose a cluster to display', value = 0, min = min(p322clusters), max = max(p322clusters), step = 1)
        updateSelectInput(session, inputId = 'p33id', label = 'Choose a sample to display', choices = names(results[['p33plots']]), selected = names(results[['p33plots']])[1])
        p34samples <- unique(gsub(".*\\|(.*)","\\1",names(results[['p34plots']])))
        p342clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p34plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p34id', label = 'Choose a sample to display', choices = p34samples, selected = p34samples[1])
        updateSliderInput(session, inputId ="p342id", label = 'Choose a cluster to display', value = 0, min = min(p342clusters), max = max(p342clusters), step = 1)
        p35samples <- unique(gsub(".*\\|(.*)","\\1",names(results[['p35plots']])))
        p352clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p35plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p35id', label = 'Choose a sample to display', choices = p35samples, selected = p35samples[1])
        updateSliderInput(session, inputId ="p352id", label = 'Choose a cluster to display', value = 0, min = min(p352clusters), max = max(p352clusters), step = 1)
        updateSelectInput(session, inputId = 'p36id', label = 'Choose a sample to display', choices = names(results[['p36plots']]), selected = names(results[['p36plots']])[1])
        updateSelectInput(session, inputId = 'p37id', label = 'Choose a sample to display', choices = names(results[['p37plots']]), selected = names(results[['p37plots']])[1])
        updateSelectInput(session, inputId = 'p373id', label = 'Choose a sample to display', choices = names(results[['p373plots']]), selected = names(results[['p373plots']])[1])
        updateSelectInput(session, inputId = 'p38id', label = 'Choose a sample to display', choices = names(results[['p38plots']]), selected = names(results[['p38plots']])[1])
        p54clusters <- unique(sort(as.numeric(as.character(gsub("CLUSTER_(.*)","\\1",names(results[['p54plots']]), ignore.case = T)))))
        updateSliderInput(session, inputId ="p54id", label = 'Choose a cluster to display', value = 0, min = min(p54clusters), max = max(p54clusters), step = 1)
        p70clusters <- sort(as.numeric(as.character(gsub("CLUSTER_(.*)","\\1",names(results[['p70plots']]), ignore.case = T))))
        updateSliderInput(session, inputId ="p70id", label = 'Choose a cluster to display', value = 0, min = min(p70clusters), max = max(p70clusters), step = 1)
        # p71clusters <- sort(as.numeric(as.character(gsub("CLUSTER_(.*)","\\1",names(results[['p71plots']]), sep = ""))))
        # updateSliderInput(session, inputId ="p71id", label = 'Choose a cluster to display', value = 0, min = min(p71clusters), max = max(p71clusters), step = 1)
        updateSelectInput(session, inputId = 'p71id', label = 'Choose a top gene to display', choices = names(results[['p71plots']]), selected = names(results[['p71plots']])[1])
        updateSelectInput(session, inputId = 'p72id', label = 'Choose a sample to display', choices = names(results[['p72plots']]), selected = names(results[['p72plots']])[1])
        updateSelectInput(session, inputId = 'p73id', label = 'Choose a sample to display', choices = names(results[['p73plots']]), selected = names(results[['p73plots']])[1])
        updateSelectInput(session, inputId = 'p74id', label = 'Choose a sample to display', choices = names(results[['p74plots']]), selected = names(results[['p74plots']])[1])
        updateSelectInput(session, inputId = 'p75id', label = 'Choose a sample to display', choices = names(results[['p75plots']]), selected = names(results[['p75plots']])[1])
        updateSelectInput(session, inputId = 'p76id', label = 'Choose a sample to display', choices = names(results[['p76plots']]), selected = names(results[['p76plots']])[1])
        p77samples <- unique(gsub(".*\\|(.*)","\\1",names(results[['p77plots']])))
        p77clusters <- sort(as.numeric(as.character(unique(gsub("CLUSTER_(.*)\\|.*","\\1",names(results[['p77plots']]))))), decreasing = F)
        updateSelectInput(session, inputId = 'p77id', label = 'Choose a sample to display', choices = p77samples, selected = p77samples[1])
        updateSliderInput(session, inputId ="p772id", label = 'Choose a cluster to display', value = 0, min = min(p77clusters), max = max(p77clusters), step = 1)
        p78samples <- unique(gsub("(.*?)\\|.*\\|.*","\\1",names(results[['p78plots']])))
        p78types <- unique(gsub(".*?\\|(.*)\\|.*","\\1",names(results[['p78plots']])))
        p78para <- unique(gsub(".*?\\|.*\\|(.*)","\\1",names(results[['p78plots']])))
        updateSelectInput(session, inputId = 'p78id', label = 'Choose a sample to display', choices = p78samples, selected = p78samples[1])
        updateSelectInput(session, inputId = 'p782id', label = 'Choose a group to display', choices = p78types, selected = p78types[1])
        updateSelectInput(session, inputId = 'p783id', label = 'Choose a parameter to display', choices = p78para, selected = p78para[1])
        p79samples <- unique(gsub("(.*)\\|.*","\\1",names(results[['p79plots']])))
        p792celltypes <- unique(gsub(".*\\|(.*)","\\1",names(results[['p79plots']])))
        updateSelectInput(session, inputId = 'p79id', label = 'Choose a sample to display', choices = p79samples, selected = p79samples[1])
        updateSelectInput(session, inputId = 'p792id', label = 'Choose a cell type to display', choices = p792celltypes, selected = p792celltypes[1])
        updateSelectInput(session, inputId = 'p80id', label = 'Choose a sample to display', choices = names(results[['p80plots']]), selected = names(results[['p80plots']])[1])
        updateSelectInput(session, inputId = 'p81id', label = 'Choose a sample to display', choices = names(results[['p81plots']]), selected = names(results[['p81plots']])[1])
        updateSelectInput(session, inputId = 'p82id', label = 'Choose a sample to display', choices = names(results[['p82plots']]), selected = names(results[['p82plots']])[1])
        updateSelectInput(session, inputId = 'p83id', label = 'Choose a sample to display', choices = names(results[['p83plots']]), selected = names(results[['p83plots']])[1])
        p84types <- unique(gsub(".*?\\|(.*)","\\1",names(results[['p84plots']])))
        p84para <- unique(gsub("(.*?)\\|.*","\\1",names(results[['p84plots']])))
        updateSelectInput(session, inputId = 'p84id', label = 'Choose a group to display', choices = p84types, selected = p84types[1])
        updateSelectInput(session, inputId = 'p842id', label = 'Choose a variable to display', choices = p84para, selected = p84para[1])
        updateSliderInput(session, inputId ="p90id", label = 'Choose a cluster to display', value = 0, min = min(as.numeric(names(results[['p90plots']]))), max = max(as.numeric(names(results[['p90plots']]))), step = 1)
        print("Completed!")
        return(results)
    })
    
    output$p1plot <- renderPlot({
        showModal(modalDialog("Plotting figure 1..", footer=NULL))
        p <- data()[['p1plots']][[input$p1id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p1plot, {
        screenshot(id="p1plot", filename = paste("1SCA_TABLE_PREFILTERING_RNA_INFO_",input$p1id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p2plot <- renderPlot({
        showModal(modalDialog("Plotting figure 2..", footer=NULL))
        p <- data()[['p2plots']][[input$p2id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p2plot, {
        screenshot(id="p2plot", filename = paste("2SCA_PLOT_PREFILTERING_RNA_INFO_VIOLIN_",input$p2id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p3plot <- renderPlot({
        if(length(grep("ggplot", class(data()[['p3plots']][[input$p3id]]), ignore.case = T)) > 0){
            showModal(modalDialog("Plotting figure 3..", footer=NULL))
            p <- data()[['p3plots']][[input$p3id]]
            removeModal()
        }else{
            text <- paste("\n   Insufficient cell number\n",
                          "     provided .")
            p <- ggplot() +
                annotate("text", x = 4, y = 25, size=8, label = text) +
                theme_void()
        }
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p3plot, {
        screenshot(id="p3plot", filename = paste("3SCA_PLOT_CELL_CYCLE_PHASES_PCA_BEFORE_CC_REGRESSION_",input$p3id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p4plot <- renderPlot({
        if(length(grep("ggplot", class(data()[['p4plots']][[input$p4id]]), ignore.case = T)) > 0){
            showModal(modalDialog("Plotting figure 4..", footer=NULL))
            p <- data()[['p4plots']][[input$p4id]]
            removeModal()
        }else{
            text <- paste("\n   Insufficient cell number\n",
                          "     provided.")
            p <- ggplot() +
                annotate("text", x = 4, y = 25, size=8, label = text) +
                theme_void()
        }
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p4plot, {
        screenshot(id="p4plot", filename = paste("4SCA_PLOT_CELL_CYCLE_PHASES_UMAP_BEFORE_CC_REGRESSION_",input$p4id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p5plot <- renderPlot({
      showModal(modalDialog("Plotting figure 5..", footer=NULL))
        if(length(grep("ggplot", class(data()[['p5plots']][[input$p5id]]), ignore.case = T)) > 0){
          p <- data()[['p5plots']][[input$p5id]]
          removeModal()
        }else{
          text <- paste("\n   Insufficient cell number\n",
                        "     provided.")
          p <- ggplot() +
            annotate("text", x = 4, y = 25, size=8, label = text) +
            theme_void()
        }
        removeModal()
        return(p)
    }, height = 800, width = 900)
    
    observeEvent(input$p5plot, {
        screenshot(id="p5plot", filename = paste("5SCA_PLOT_METRICS_UMAP_BEFORE_CC_REGRESSION_",input$p5id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p6plot <- renderPlot({
        if(length(grep("ggplot", class(data()[['p6plots']][[input$p6id]]), ignore.case = T)) > 0){
            showModal(modalDialog("Plotting figure 6..", footer=NULL))
            p <- data()[['p6plots']][[input$p6id]]
            removeModal()
        }else{
            text <- paste("\n   No cell cycle regression\n",
                          "     has been carried out.")
            p <- ggplot() +
                annotate("text", x = 4, y = 25, size=8, label = text) +
                theme_void()
        }
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p6plot, {
        screenshot(id="p6plot", filename = paste("6SCA_PLOT_CELL_CYCLE_PHASES_PCA_POST_CC_REGRESSION_",input$p6id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p7plot <- renderPlot({
        if(length(grep("ggplot", class(data()[['p7plots']][[input$p7id]]), ignore.case = T)) > 0){
            showModal(modalDialog("Plotting figure 7..", footer=NULL))
            p <- data()[['p7plots']][[input$p7id]]
            removeModal()
        }else{
            text <- paste("\n   No cell cycle regression\n",
                          "     has been carried out.")
            p <- ggplot() +
                annotate("text", x = 4, y = 25, size=8, label = text) +
                theme_void()
        }
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p7plot, {
        screenshot(id="p7plot", filename = paste("7SCA_PLOT_CELL_CYCLE_PHASES_UMAP_POST_CC_REGRESSION_",input$p7id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p8plot <- renderPlot({
        showModal(modalDialog("Plotting figure 8..", footer=NULL))
        p <- data()[['p8plots']][[input$p8id]]
        removeModal()
        return(p)
    }, height = 800, width = 900)
    
    observeEvent(input$p8plot, {
        screenshot(id="p8plot", filename = paste("8SCA_PLOT_METRICS_UMAP_POST_CC_REGRESSION_",input$p8id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p9plot <- renderPlot({
        showModal(modalDialog("Plotting figure 9..", footer=NULL))
        p <- data()[['p9plots']][[input$p9id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p9plot, {
        screenshot(id="p9plot", filename = paste("9SCA_PLOT_POSTFILTERING_RNA_INFO_VIOLIN_",input$p9id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p11plot <- renderPlot({
        showModal(modalDialog("Plotting figure 11..", footer=NULL))
        p <- data()[['p11plots']][[input$p11id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p11plot, {
        screenshot(id="p11plot", filename = paste("11SCA_PLOT_POSTNORMALIZED_VARIABLE_FEATURE_",input$p11id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p12plot <- renderPlot({
        showModal(modalDialog("Plotting figure 12..", footer=NULL))
        p <- data()[['p12plots']][[input$p12id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p12plot, {
        screenshot(id="p12plot", filename = paste("12SCA_PLOT_TOP_FEATURES_ASSOC_PC12_",input$p12id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p13plot <- renderPlot({
        showModal(modalDialog("Plotting figure 13..", footer=NULL))
        p <- data()[['p13plots']][[input$p13id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p13plot, {
        screenshot(id="p13plot", filename = paste("13SCA_PLOT_VISUALIZE_PC12_",input$p13id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p14plot <- renderPlot({
        showModal(modalDialog("Plotting figure 14..", footer=NULL))
        p <- data()[['p14plots']][[input$p14id]]
        removeModal()
        return(p)
    }, height = 1000, width = 900)
    
    observeEvent(input$p14plot, {
        screenshot(id="p14plot", filename = paste("14SCA_PLOT_TOP_FEATURES_ASSOC_PC15_",input$p14id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p15plot <- renderPlot({
        showModal(modalDialog("Plotting figure 15..", footer=NULL))
        p <- data()[['p15plots']][[input$p15id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p15plot, {
        screenshot(id="p15plot", filename = paste("15SCA_PLOT_UMAP_AUTOCLUSTER_",input$p15id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p16plot <- renderPlot({
        showModal(modalDialog("Plotting figure 16..", footer=NULL))
        p <- data()[['p16plots']][[input$p16id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p16plot, {
        screenshot(id="p16plot", filename = paste("16SCA_PLOT_TSNE_AUTOCLUSTER_",input$p16id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p17table <- renderDT(data()[['p17data']][[input$p17id]],
                                filter = "top",
                                style="bootstrap",
                                rownames = F,
                                options = list(pageLength = 10))
    
    output$p18plot <- renderPlot({
        showModal(modalDialog("Plotting figure 18..", footer=NULL))
        p <- data()[['p18plots']][[input$p18id]]
        removeModal()
        return(p)
    }, height = 800, width = 900)
    
    observeEvent(input$p18plot, {
        screenshot(id="p18plot", filename = paste("18SCA_PLOT_HEATMAP_MEDIAN_EXPR_CLUSTER_",input$p18id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p19plot <- renderPlot({
        showModal(modalDialog("Plotting figure 19..", footer=NULL))
        p1 <- data()[['p19plots']][[input$p19id]]
        p2 <- data()[['p192plots']][[input$p19id]]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 450, width = 900)
    
    observeEvent(input$p19plot, {
        screenshot(id="p19plot", filename = paste("19SCA_PLOT_AVE_LOGFC_TOP_GENES_TSNE_",input$p19id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p193plot <- renderPlot({
        showModal(modalDialog("Plotting figure 20..", footer=NULL))
        p1 <- data()[['p19plots']][[input$p193id]]
        p2 <- data()[['p193plots']][[input$p193id]]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 450, width = 900)
    
    observeEvent(input$p193plot, {
        screenshot(id="p193plot", filename = paste("192SCA_scRNASEQ_AVE_LOGFC_TOP_GENES_UMAP_",input$p193id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p20plot <- renderPlot({
        showModal(modalDialog("Plotting figure 20..", footer=NULL))
        p <- data()[['p20plots']][[input$p20id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p20plot, {
        screenshot(id="p20plot", filename = paste("20SCA_PLOT_RIDGE_TOP_FC_SIG_GENES_IN_CLUSTERS_",input$p20id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p21plot <- renderPlot({
        showModal(modalDialog("Plotting figure 21..", footer=NULL))
        p <- data()[['p21plots']][[input$p21id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p21plot, {
        screenshot(id="p21plot", filename = paste("21SCA_PLOT_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",input$p21id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p22plot <- renderPlot({
        showModal(modalDialog("Plotting figure 22..", footer=NULL))
        p <- data()[['p22plots']][[paste("CLUSTER_",input$p222id,"|",input$p22id, sep = "")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p22plot, {
        screenshot(id="p22plot", filename = paste("22SCA_PLOT_VIOLIN_TOP_MARKERS_IN_CLUSTERS_",input$p22id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p23plot <- renderPlot({
        showModal(modalDialog("Plotting figure 23..", footer=NULL))
        p <- data()[['p23plots']][[input$p23id]]
        removeModal()
        return(p)
    }, height = 900, width = 900)
    
    observeEvent(input$p23plot, {
        screenshot(id="p23plot", filename = paste("23SCA_PLOT_HEATMAP_TOP_MARKERS_IN_CLUSTERS_",input$p23id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p24plot <- renderPlot({
        showModal(modalDialog("Plotting figure 24..", footer=NULL))
        p <- data()[['p24plots']][[input$p24id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p24plot, {
        screenshot(id="p24plot", filename = paste("24SCA_PLOT_UMAP_AUTO_CELL_TYPE_IDENTIFICATION_",input$p24id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p25plot <- renderPlot({
        showModal(modalDialog("Plotting figure 25..", footer=NULL))
        p <- data()[['p25plots']][[input$p25id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p25plot, {
        screenshot(id="p25plot", filename = paste("25SCA_PLOT_TSNE_AUTO_CELL_TYPE_IDENTIFICATION_",input$p25id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p26plot <- renderPlot({
        showModal(modalDialog("Plotting figure 26..", footer=NULL))
        p <- data()[['p26plots']][[paste("CLUSTER_",input$p262id,"|",input$p26id, sep = "")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p26plot, {
        screenshot(id="p26plot", filename = paste("26SCA_PLOT_BARPLOT_TOP_PATHWAYS_FROM_TOPFC200GENES_CLUSTER_",input$p26id,"_",input$p262id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p27table <- renderDT(data()[['p27data']][[input$p27id]],
                                filter = "top",
                                style="bootstrap",
                                rownames = F,
                                options = list(pageLength = 10, scrollX = T))
    
    output$p28plot <- renderPlot({
        showModal(modalDialog("Plotting figure 28..", footer=NULL))
        p <- data()[['p28plots']][[paste("CLUSTER_",input$p282id,"|",input$p28id, sep = "")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p28plot, {
        screenshot(id="p28plot", filename = paste("28SCA_PLOT_DOTPLOT_TOP_PATHWAYS_FROM_TOPFC200GENES_CLUSTER_",input$p28id,"_",input$p282id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p29plot <- renderPlot({
        showModal(modalDialog("Plotting figure 29..", footer=NULL))
        p <- data()[['p29plots']][[paste("CLUSTER_",input$p292id,"|",input$p29id, sep = "")]]
        removeModal()
        return(p)
    }, height = 300, width = 900)
    
    observeEvent(input$p29plot, {
        screenshot(id="p29plot", filename = paste("29SCA_PLOT_GENE_CONCEPT_NETWORK_CLUSTER_",input$p29id,"_",input$p292id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p30plot <- renderPlot({
        showModal(modalDialog("Plotting figure 30..", footer=NULL))
        p <- data()[['p30plots']][[paste("CLUSTER_",input$p302id,"|",input$p30id, sep = "")]]
        removeModal()
        return(p)
    }, height = 300, width = 900)
    
    observeEvent(input$p30plot, {
        screenshot(id="p30plot", filename = paste("30SCA_PLOT_NETWORK_ENRICHMENT_TERMS_",input$p30id,"_",input$p302id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p31plot <- renderPlot({
        if(!is.null(data()[['p31plots']][[paste("CLUSTER_",input$p312id,"|",input$p31id, sep = "")]])){
        showModal(modalDialog("Plotting figure 31..", footer=NULL))
        p <- data()[['p31plots']][[paste("CLUSTER_",input$p312id,"|",input$p31id, sep = "")]]
        removeModal()
    }else{
        text <- paste("\n   No information is\n",
                      "     found for this cluster number.")
        p <- ggplot() +
            annotate("text", x = 4, y = 25, size=8, label = text) +
            theme_void()
    }
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p31plot, {
        screenshot(id="p31plot", filename = paste("31SCA_PLOT_HEATMAP_NETWORK_TERMS_GENES_CLUSTER_",input$p31id,"_",input$p312id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p32plot <- renderPlot({
        showModal(modalDialog("Plotting figure 32..", footer=NULL))
        p <- data()[['p32plots']][[paste("CLUSTER_",input$p322id,"|",input$p32id, sep = "")]]
        removeModal()
        return(p)
    }, height = 1000, width = 900)
    
    observeEvent(input$p32plot, {
        screenshot(id="p32plot", filename = paste("32SCA_PLOT_ENRICHMENT_PROPORTIONS_GENES_CLUSTER_",input$p32id,"_",input$p322id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p33plot <- renderPlot({
        showModal(modalDialog("Plotting figure 33..", footer=NULL))
        p <- data()[['p33plots']][[input$p33id]]
        removeModal()
        return(p)
    }, height = 1000, width = 900)
    
    observeEvent(input$p33plot, {
        screenshot(id="p33plot", filename = paste("33SCA_PLOT_KEGG_ENRICHMENT_MAP_PROPORTION_CLUSTER_",input$p33id,"_",input$p332id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p34plot <- renderPlot({
        showModal(modalDialog("Plotting figure 34..", footer=NULL))
        p <- data()[['p34plots']][[paste("CLUSTER_",input$p342id,"|",input$p34id, sep = "")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p34plot, {
        screenshot(id="p34plot", filename = paste("34SCA_PLOT_UPSET_PLOT_OVERLAPS_BETWEEN_GENE_SETS_",input$p34id,"_",input$p342id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p35plot <- renderPlot({
        showModal(modalDialog("Plotting figure 35..", footer=NULL))
        p <- data()[['p35plots']][[paste("CLUSTER_",input$p352id,"|",input$p35id, sep = "")]]
        removeModal()
        return(p)
    }, height = 500, width = 900)
    
    observeEvent(input$p35plot, {
        screenshot(id="p35plot", filename = paste("35SCA_PLOT_GSEA_GENE_ENRICHMENT_SCORE_",input$p35id,"_",input$p352id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p36plot <- renderPlot({
        showModal(modalDialog("Plotting figure 36..", footer=NULL))
        p <- data()[['p36plots']][[input$p36id]]
        removeModal()
        return(p)
    }, height = 1200, width = 900)
    
    observeEvent(input$p36plot, {
        screenshot(id="p36plot", filename = paste("36SCA_PLOT_HEATMAP_MEDIAN_EXPR_CELL_TYPE_TOP_FOLDCHANGE_SIG_GENES_",input$p36id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p37plot <- renderPlot({
        showModal(modalDialog("Plotting figure 37..", footer=NULL))
        p1 <- data()[['p37plots']][[input$p37id]]
        p2 <- data()[['p372plots']][[input$p37id]]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 450, width = 900)
    
    observeEvent(input$p37plot, {
        screenshot(id="p37plot", filename = paste("37SCA_PLOT_AVE_LOGFC_TOP_GENES_TSNE_",input$p37id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p373plot <- renderPlot({
      showModal(modalDialog("Plotting figure 373..", footer=NULL))
      p1 <- data()[['p37plots']][[input$p373id]]
      p2 <- data()[['p373plots']][[input$p373id]]
      removeModal()
      print(grid.arrange(p1, p2, ncol = 2))
    }, height = 450, width = 900)
    
    observeEvent(input$p373plot, {
      screenshot(id="p373plot", filename = paste("38SCA_PLOT_AVE_LOGFC_TOP_GENES_UMAP_",input$p373id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p38plot <- renderPlot({
        showModal(modalDialog("Plotting figure 38..", footer=NULL))
        p <- data()[['p38plots']][[input$p38id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p38plot, {
        screenshot(id="p38plot", filename = paste("38SCA_PLOT_AVE_LOGFC_TOP_GENES_UMAP_",input$p38id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p39plot <- renderPlot({
        showModal(modalDialog("Plotting figure 39..", footer=NULL))
        p <- data()[['p39plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p39plot, {
        screenshot(id="p39plot", filename = paste("39SCA_PLOT_UMAP_BEFORE_AFTER_SAMPLE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p40plot <- renderPlot({
        showModal(modalDialog("Plotting figure 40..", footer=NULL))
        p <- data()[['p40plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p40plot, {
        screenshot(id="p40plot", filename = paste("40SCA_PLOT_TSNE_BEFORE_AFTER_SAMPLE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p41plot <- renderPlot({
        showModal(modalDialog("Plotting figure 41..", footer=NULL))
        p <- data()[['p41plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p41plot, {
        screenshot(id="p41plot", filename = paste("41SCA_PLOT_PCA_BEFORE_AFTER_SAMPLE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p42plot <- renderPlot({
        showModal(modalDialog("Plotting figure 42..", footer=NULL))
        p <- data()[['p42plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p42plot, {
        screenshot(id="p42plot", filename = paste("42SCA_PLOT_UMAP_BEFORE_AFTER_GROUP_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p43plot <- renderPlot({
        showModal(modalDialog("Plotting figure 43..", footer=NULL))
        p <- data()[['p43plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p43plot, {
        screenshot(id="p43plot", filename = paste("43SCA_PLOT_TSNE_BEFORE_AFTER_GROUP_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p44plot <- renderPlot({
        showModal(modalDialog("Plotting figure 44..", footer=NULL))
        p <- data()[['p44plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p44plot, {
        screenshot(id="p44plot", filename = paste("44SCA_PLOT_PCA_BEFORE_AFTER_GROUP_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p45plot <- renderPlot({
        showModal(modalDialog("Plotting figure 45..", footer=NULL))
        p <- data()[['p45plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p45plot, {
        screenshot(id="p45plot", filename = paste("45SCA_PLOT_UMAP_BEFORE_AFTER_BATCH_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p46plot <- renderPlot({
        showModal(modalDialog("Plotting figure 46..", footer=NULL))
        p <- data()[['p46plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p46plot, {
        screenshot(id="p46plot", filename = paste("46SCA_PLOT_TSNE_BEFORE_AFTER_BATCH_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p47plot <- renderPlot({
        showModal(modalDialog("Plotting figure 47..", footer=NULL))
        p <- data()[['p47plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p47plot, {
        screenshot(id="p47plot", filename = paste("47SCA_PLOT_PCA_BEFORE_AFTER_BATCH_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p48plot <- renderPlot({
        showModal(modalDialog("Plotting figure 48..", footer=NULL))
        p <- data()[['p48plots']]
        removeModal()
        return(p)
    }, height = 600, width = 1400)
    
    observeEvent(input$p48plot, {
        screenshot(id="p48plot", filename = paste("48SCA_PLOT_TOP_FEATURES_ASSOC_PC12_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p49plot <- renderPlot({
        showModal(modalDialog("Plotting figure 49..", footer=NULL))
        p <- data()[['p49plots']]
        removeModal()
        return(p)
    }, height = 1600, width = 1200)
    
    observeEvent(input$p49plot, {
        screenshot(id="p49plot", filename = paste("49SCA_PLOT_HEATMAP_PC1_15_TOP_FEATURES_INTEGRATED_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p50plot <- renderPlot({
        showModal(modalDialog("Plotting figure 50..", footer=NULL))
        p <- data()[['p50plots']]
        removeModal()
        return(p)
    }, height = 700, width = 1200)
    
    observeEvent(input$p50plot, {
        screenshot(id="p50plot", filename = paste("50SCA_PLOT_UMAP_AUTOCLUSTER_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p51plot <- renderPlot({
        showModal(modalDialog("Plotting figure 51..", footer=NULL))
        p <- data()[['p51plots']]
        removeModal()
        return(p)
    }, height = 700, width = 1200)
    
    observeEvent(input$p51plot, {
        screenshot(id="p51plot", filename = paste("51SCA_PLOT_TSNE_AUTOCLUSTER_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p52table <- renderDT(data()[['p52data']][[input$p52id]],
                                filter = "top",
                                style="bootstrap",
                                rownames = F,
                                options = list(pageLength = 10))
    
    output$p53plot <- renderPlot({
        showModal(modalDialog("Plotting figure 53..", footer=NULL))
        p <- data()[['p53plots']]
        removeModal()
        return(p)
    }, height = 1400, width = 1200)
    
    observeEvent(input$p53plot, {
        screenshot(id="p53plot", filename = paste("53SCA_PLOT_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_","_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p54plot <- renderPlot({
        showModal(modalDialog("Plotting figure 54..", footer=NULL))
        p <- data()[['p54plots']][[as.character(input$p54id)]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p54plot, {
        screenshot(id="p54plot", filename = paste("54SCA_PLOT_VIOLIN_TOP_GENES_CLUSTER_",input$p54id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p55plot <- renderPlot({
        showModal(modalDialog("Plotting figure 55..", footer=NULL))
        p <- data()[['p55plots']]
        removeModal()
        return(p)
    }, height = 1400, width = 1200)
    
    observeEvent(input$p55plot, {
        screenshot(id="p55plot", filename = paste("55SCA_PLOT_HEATMAP_TOP_MARKERS_IN_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p56plot <- renderPlot({
        showModal(modalDialog("Plotting figure 56..", footer=NULL))
        p <- data()[['p56plots']]
        removeModal()
        return(p)
    }, height = 1200, width = 1400)
    
    observeEvent(input$p56plot, {
        screenshot(id="p56plot", filename = paste("56SCA_PLOT_UMAP_AUTO_CELL_TYPE_IDENTIFICATION_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p57plot <- renderPlot({
        showModal(modalDialog("Plotting figure 57..", footer=NULL))
        p <- data()[['p57plots']]
        removeModal()
        return(p)
    }, height = 1200, width = 1400)
    
    observeEvent(input$p57plot, {
        screenshot(id="p57plot", filename = paste("57SCA_PLOT_TSNE_AUTO_CELL_TYPE_IDENTIFICATION_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p58plot <- renderPlot({
        showModal(modalDialog("Plotting figure 58..", footer=NULL))
        p <- data()[['p58plots']]
        removeModal()
        return(p)
    }, height = 500, width = 900)
    
    observeEvent(input$p58plot, {
        screenshot(id="p58plot", filename = paste("59SCA_PLOT_UMAP_BY_GROUP_AUTO_CELL_TYPE_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p59plot <- renderPlot({
        showModal(modalDialog("Plotting figure 59..", footer=NULL))
        p <- data()[['p59plots']]
        removeModal()
        return(p)
    }, height = 1600, width = 1200)
    
    observeEvent(input$p59plot, {
        screenshot(id="p59plot", filename = paste("59SCA_PLOT_HEATMAP_MEDIAN_EXPR_CLUSTER_TOP_FOLDCHANGE_SIG_GENES_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p60plot <- renderPlot({
        showModal(modalDialog("Plotting figure 60..", footer=NULL))
      p1 <- data()[['p60plots']]
      p2 <- data()[['p602plots']]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 600, width = 1400)
    
    observeEvent(input$p60plot, {
        screenshot(id="p60plot", filename = paste("60SCA_PLOT_AVE_LOGFC_TOP_GENES_TSNE_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p602plot <- renderPlot({
        showModal(modalDialog("Plotting figure 61..", footer=NULL))
        p1 <- data()[['p60plots']]
        p2 <- data()[['p603plots']]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 600, width = 1400)
    
    observeEvent(input$p602plot, {
        screenshot(id="p602plot", filename = paste("602SCA_PLOT_AVE_LOGFC_TOP_GENES_UMAP_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p61plot <- renderPlot({
        showModal(modalDialog("Plotting figure 61..", footer=NULL))
        p <- data()[['p61plots']]
        removeModal()
        return(p)
    }, height = 1000, width = 1400)
    
    observeEvent(input$p61plot, {
        screenshot(id="p61plot", filename = paste("61SCA_PLOT_RIDGE_TOP_FC_SIG_GENES_IN_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p62plot <- renderPlot({
        showModal(modalDialog("Plotting figure 62..", footer=NULL))
        p <- data()[['p62plots']]
        removeModal()
        return(p)
    }, height = 1200, width = 1200)
    
    observeEvent(input$p62plot, {
        screenshot(id="p62plot", filename = paste("62SCA_PLOT_NODE_SUMMARY_SAMPLE_PROPORTION_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p63plot <- renderPlot({
        showModal(modalDialog("Plotting figure 63..", footer=NULL))
        p <- data()[['p63plots']]
        removeModal()
        return(p)
    }, height = 1000, width = 1400)
    
    observeEvent(input$p63plot, {
        screenshot(id="p63plot", filename = paste("63SCA_PLOT_CELL_TYPE_SUMMARY_SAMPLE_PROPORTION_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p64plot <- renderPlot({
        showModal(modalDialog("Plotting figure 64..", footer=NULL))
        p <- data()[['p64plots']]
        removeModal()
        return(p)
    }, height = 1400, width = 1200)
    
    observeEvent(input$p64plot, {
        screenshot(id="p64plot", filename = paste("64SCA_PLOT_HEATMAP_MEDIAN_EXPR_CELL_TYPE_TOP_FOLDCHANGE_SIG_GENES_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p65plot <- renderPlot({
        showModal(modalDialog("Plotting figure 65..", footer=NULL))
        p1 <- data()[['p65plots']]
        p2 <- data()[['p652plots']]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 700, width = 1400)
    
    observeEvent(input$p65plot, {
        screenshot(id="p65plot", filename = paste("65SCA_PLOT_AVE_LOGFC_TOP_GENES_TSNE_CELL_TYPES_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p652plot <- renderPlot({
        showModal(modalDialog("Plotting figure 652..", footer=NULL))
        p1 <- data()[['p65plots']]
        p2 <- data()[['p653plots']]
        removeModal()
        print(grid.arrange(p1, p2, ncol = 2))
    }, height = 700, width = 1400)
    
    observeEvent(input$p652plot, {
        screenshot(id="p652plot", filename = paste("65SCA_PLOT_AVE_LOGFC_TOP_GENES_UMAP_CELL_TYPES_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p66plot <- renderPlot({
        showModal(modalDialog("Plotting figure 66..", footer=NULL))
        p <- data()[['p66plots']]
        removeModal()
        return(p)
    }, height = 1000, width = 1400)
    
    observeEvent(input$p66plot, {
        screenshot(id="p66plot", filename = paste("66SCA_PLOT_RIDGE_TOP_FC_SIG_GENES_CELL_TYPES_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p67plot <- renderPlot({
        showModal(modalDialog("Plotting figure 67..", footer=NULL))
        p <- data()[['p67plots']]
        removeModal()
        return(p)
    }, height = 1200, width = 1200)
    
    observeEvent(input$p67plot, {
        screenshot(id="p67plot", filename = paste("67SCA_PLOT_CELL_TYPE_SAMPLE_PROPORTION_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p68table <- renderDT(data()[['p68data']][[input$p68id]],
                                filter = "top",
                                style="bootstrap",
                                rownames = F,
                                options = list(pageLength = 10))
    
    output$p69plot <- renderPlot({
        showModal(modalDialog("Plotting figure 69..", footer=NULL))
        p <- data()[['p69plots']]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p69plot, {
        screenshot(id="p69plot", filename = paste("69SCA_PLOT_CLUSTERS_BY_GROUP_TOP1_FEATURE_DIFFERENCE_COLOR_EXPR_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p70plot <- renderPlot({
        showModal(modalDialog("Plotting figure 70..", footer=NULL))
        p <- data()[['p70plots']][[paste("CLUSTER_",input$p70id,sep = "")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p70plot, {
        screenshot(id="p70plot", filename = paste("70SCA_PLOT_VIOLIN_TOP_MARKERS_DE_GROUPS_CLUSTER_",input$p70id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p71plot <- renderPlot({
        showModal(modalDialog("Plotting figure 71..", footer=NULL))
        p <- data()[['p71plots']][[input$p71id]]
        removeModal()
        return(p)
    }, height = 450, width = 900)
    
    observeEvent(input$p71plot, {
        screenshot(id="p71plot", filename = paste("71SCA_PLOT_UMAP_DENSITY_TOP_1_MARKER_BY_DE_GROUP_CLUSTER_",gsub("\\s+","",input$p71id),"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p72plot <- renderPlot({
        showModal(modalDialog("Plotting figure 72..", footer=NULL))
        p <- data()[['p72plots']][[input$p72id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p72plot, {
        screenshot(id="p72plot", filename = paste("72SCA_PLOT_UMAP_TRAJECTORY_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p73plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 73..", footer=NULL))
        p <- data()[['p73plots']][[input$p73id]]
        removeModal()
        return(p)
    })
    
    output$p74plot <- renderPlot({
        showModal(modalDialog("Plotting figure 74..", footer=NULL))
        p <- data()[['p74plots']][[input$p74id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p74plot, {
        screenshot(id="p74plot", filename = paste("74SCA_PLOT_UMAP_TRAJECTORY_CELL_TYPE_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p75plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 75..", footer=NULL))
        p <- data()[['p75plots']][[input$p75id]]
        removeModal()
        return(p)
    })
    
    output$p76plot <- renderPlot({
        showModal(modalDialog("Plotting figure 76..", footer=NULL))
        p <- data()[['p76plots']][[input$p76id]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p76plot, {
        screenshot(id="p76plot", filename = paste("76SCA_PLOT_UMAP_TRAJECTORY_PARTITION_",input$p76id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p77plot <- renderPlot({
        showModal(modalDialog("Plotting figure 77..", footer=NULL))
        p <- data()[['p77plots']][[paste("CLUSTER_",input$p772id,"|",input$p77id,sep = "")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p77plot, {
        screenshot(id="p77plot", filename = paste("77SCA_PLOT_UMAP_PSEUDOTIME_ROOT_CLUSTER_",input$p77id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p78plot <- renderPlot({
        showModal(modalDialog("Plotting figure 78..", footer=NULL))
        p <- data()[['p78plots']][[paste(input$p78id,input$p782id,input$p783id,sep = "|")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p78plot, {
        screenshot(id="p78plot", filename = paste("78SCA_PLOT_",input$p782id,"_VS_",input$p783id,"_",input$p78id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p79plot <- renderPlot({
        showModal(modalDialog("Plotting figure 79..", footer=NULL))
        p <- data()[['p79plots']][[paste(input$p79id,input$p792id,sep = "|")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p79plot, {
        screenshot(id="p79plot", filename = paste("79SCA_PLOT_",input$p792id,"_",input$p79id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p80plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 80..", footer=NULL))
        p <- data()[['p80plots']][[input$p80id]]
        removeModal()
        return(p)
    })
    
    output$p81plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 81..", footer=NULL))
        p <- data()[['p81plots']][[input$p81id]]
        removeModal()
        return(p)
    })
    
    output$p82plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 82..", footer=NULL))
        p <- data()[['p82plots']][[input$p82id]]
        removeModal()
        return(p)
    })
    
    output$p83plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 83..", footer=NULL))
        p <- data()[['p83plots']][[input$p83id]]
        removeModal()
        return(p)
    })
    
    output$p84plot <- renderPlot({
        showModal(modalDialog("Plotting figure 84..", footer=NULL))
        p <- data()[['p84plots']][[paste(input$p842id,input$p84id,sep = "|")]]
        removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p84plot, {
        screenshot(id="p84plot", filename = paste("84SCA_PLOT_CELL_TRAJECTORY_BY_TIMEPOINT_STATE_",input$p842id,"_VS_",input$p84id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p85plot <- renderPlot({
        showModal(modalDialog("Plotting figure 85..", footer=NULL))
        p <- data()[['p85plots']]
        removeModal()
        return(p)
    }, height = 900, width = 1400)
    
    observeEvent(input$p85plot, {
        screenshot(id="p85plot", filename = paste("85SCA_PLOT_UMAP_TRAJECTORY_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p86plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 86..", footer=NULL))
        p <- data()[['p86plots']]
        removeModal()
        return(p)
    })
    
    output$p87plot <- renderPlot({
        showModal(modalDialog("Plotting figure 87..", footer=NULL))
        p <- data()[['p87plots']]
        removeModal()
        return(p)
    }, height = 900, width = 1400)
    
    observeEvent(input$p87plot, {
        screenshot(id="p87plot", filename = paste("87SCA_PLOT_UMAP_TRAJECTORY_CELL_TYPE_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p88plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 88..", footer=NULL))
        p <- data()[['p88plots']]
        removeModal()
        return(p)
    })
    
    output$p89plot <- renderPlot({
        showModal(modalDialog("Plotting figure 89..", footer=NULL))
        p <- data()[['p89plots']]
        removeModal()
        return(p)
    }, height = 1000, width = 1400)
    
    observeEvent(input$p89plot, {
        screenshot(id="p89plot", filename = paste("89SCA_PLOT_UMAP_TRAJECTORY_PARTITION_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p90plot <- renderPlot({
        showModal(modalDialog("Plotting figure 90..", footer=NULL))
      p <- data()[['p90plots']][[as.character(input$p90id)]]
      removeModal()
        return(p)
    }, height = 600, width = 900)
    
    observeEvent(input$p90plot, {
        screenshot(id="p90plot", filename = paste("90SCA_PLOT_UMAP_PSEUDOTIME_ROOT_CLUSTER_",input$p90id,"_",data()[['project_name']], sep = ""), scale = 2)
    })
    
    output$p91plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 91..", footer=NULL))
        p <- data()[['p91plots']]
        removeModal()
        return(p)
    })
    
    output$p92plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 92..", footer=NULL))
        p <- data()[['p92plots']]
        removeModal()
        return(p)
    })
    
    output$p93plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 93..", footer=NULL))
        p <- data()[['p93plots']]
        removeModal()
        return(p)
    })
    
    output$p94plot <- renderPlotly({
        showModal(modalDialog("Plotting figure 94..", footer=NULL))
        p <- data()[['p94plots']]
        removeModal()
        return(p)
    })
    
})


