###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: CyTOF Analysis Pipeline
# Version: V1.2.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last-update date: 2023-07-20
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
library("remotes")
library("flowCore")
library("ggplot2")
library("reshape2")
library("limma")
library("gplots")
library("plot3D")
library("Biobase")
library("FlowSOM")
library("ConsensusClusterPlus")
library("Rtsne")
library("wesanderson")
library("viridis")
library("webshot")
library("Seurat")
library("plotly")
library("lattice")
library("latticeExtra")
library("wesanderson")
library("ggpubr")
library("plyr")
library("akmedoids")
library("umap")
library("patchwork")

Sys.setenv("DISPLAY"=":0")

source("DB/SCA_CyTOF_RShiny_Functions_V1.0.0.R")
color_conditions <- color_ini()

# example1 <- "DB/SCA_CyTOF_Example_From_10X.zip"
example2 <- "DB/SCA_CyTOF_Metadata_Example.csv"

ctime <- format(Sys.time(), format = "%Y%m%d%H%M%S", tz = "Europe/Stockholm")

shinyServer(function(input, output, session) {
  values <- reactiveValues(proceed = 0, proceed2 = 0, proceed3 = 0)
  
  observeEvent(input$submit & input$termscheck,{
    
    if(input$termscheck == FALSE & input$submit){
      showModal(modalDialog("Please agree to our terms and conditions before you proceed", footer=NULL, easyClose = T))
      values$proceed <- 0
    }
    
    if(input$termscheck == TRUE & input$submit){
        updateTabsetPanel(session, "nav",selected = "section1")
        values$proceed <- 1
    }
  })
  
  output$downloadExample2 <- downloadHandler(
    filename = "SCA_CyTOF_Metadata_Example.csv",
    content = function(file) {
      file.copy(example2, file)
    },
    contentType = "text/csv"
  )
  
  results <- reactiveValues()
  
  observe({
    req(input$submit)
    req(input$files)
    req(values$proceed == 1)
    showModal(modalDialog("Initialising..", footer=NULL))
    # results <- readRDS("DB/RESULTS_EXAMPLE.RDS")
    inFile <- input$files
    print(input$phenodata$datapath)
    pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "CYTOF", isDir = T)
    print(input$files)
    print(input$files$datapath)
    project_name <- ifelse(input$project_name == "","My_Project",gsub("\\s+|\\(|\\)|-|\\/|\\?","",input$project_name))
    print(project_name)
    cdir <- paste(getwd(),"/",project_name,ctime, "/",sep = "")
    dir.create(cdir)
    sample_files <- unzip(inFile$datapath, exdir = cdir)
    sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
    input_dir <- gsub("(.*\\/).*$","\\1",sample_files[1], ignore.case = T)
    print(input_dir)

    files_submit <- list.files(pattern = "*.fcs", path = input_dir, full.names = T, ignore.case = T)
    
    # removeModal()
    
    #######################################################################################################################################
    # showModal(modalDialog("Reading FCS files..", footer=NULL))
    
    data <- NULL
    pheno_data$SID <- NULL
    
    for(i in 1:nrow(pheno_data)){
      if (length(grep(pheno_data[i,grep("FILE", colnames(pheno_data),ignore.case = T)], files_submit, ignore.case = T)) > 0){
        # current_name <- toupper(gsub("(.*)\\.FCS","\\1",gsub(".*\\/","",pheno_data[i,grep("FILE", colnames(pheno_data),ignore.case = T)]),ignore.case = T))
        current_name <- paste(pheno_data[i,"SAMPLE_ID"],pheno_data[i,"GROUP"], sep = "_")
        pheno_data[i,"SID"] <- current_name
        current <- read.FCS(files_submit[grep(pheno_data[i,grep("FILE", colnames(pheno_data),ignore.case = T)], files_submit, ignore.case = T)], transformation = FALSE, truncate_max_range = FALSE)
        col_current <- as.character(parameters(current)$desc)
        if(length(which(is.na(col_current))) > 0){
          col_current[which(is.na(col_current))] <- as.character(parameters(current)$name)[which(is.na(col_current))]
        }
        current <- data.frame(FILE = current_name, exprs(current))
        colnames(current) <- c("FILE", col_current)
        data <- rbind(data, current)
      }
    }
    
    sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SID)))
    names(sample_colors) <- unique(pheno_data$SID)
    
    colnames(data) <- gsub("^X\\.","",colnames(data), perl = TRUE)
    colnames(data) <- gsub("^CR5","CXCR5",colnames(data), perl = TRUE)
    data$FILE <- gsub(".*?/","",data$FILE, perl = TRUE)
    data$FILE <- gsub("^\\s+|\\s+$|\\s+|\\s+$|\\.fcs","",data$FILE, perl = TRUE)
    data$FILE <- toupper(data$FILE)
    colnames(data) <- gsub("\\.|_|-|\\s+","_",colnames(data), perl = TRUE)
    colnames(data)[grep("FILE", colnames(data), ignore.case = T, invert = T)] <- toupper(colnames(data))[grep("FILE", colnames(data), ignore.case = T, invert = T)]
    
    markers <- NULL
    markers <- sort(unique(colnames(data)[grep("FILE|LENGTH|DNA.*[0-9]+|TIME", colnames(data), ignore.case = T, invert = T)]))
    
    results$data <- data
    results$pheno_data <- pheno_data
    results$sample_colors <- sample_colors
    results$project_name <- project_name
    results$input_dir <- input_dir
    results$allmarkers <- markers

    removeModal()
    
    updatePickerInput(session,
      inputId = "selectmarkers", label = "Select markers",
      choices = markers,
      options = list(`actions-box` = TRUE, `selected-text-format` = "count > 2",
                     `count-selected-text` = "{0}/{1} selected"))
    
})
  
  observeEvent(input$submitmarkers,{

    if(length(input$selectmarkers) == 0 & input$submitmarkers){
      showModal(modalDialog("Please select the markers to be included in the downstream analysis before you proceed", footer=NULL, easyClose = T))
    }

    if(length(input$selectmarkers) < 4 & input$submitmarkers){
      showModal(modalDialog("Please select at least three markers before you proceed", footer=NULL, easyClose = T))
    }

    if(length(input$selectmarkers) >= 4 & input$submitmarkers){
      updateTabsetPanel(session, "nav",selected = "section2")
      values$proceed2 <- 1
    }
  })
  
  observe({
    req(!is.null(input$selectmarkers))
    req(values$proceed2 == 1)
    req(input$submitmarkers)
    req(is.null(results$som))
    
    showModal(modalDialog("Conducting analysis based on selected markers..", footer=NULL))
    data <- NULL
    data <- results$data
    results$data <- NULL
    cofactor <- 5
    files <- NULL
    files <- data$FILE
    data <- data[,grep("FILE|time|cell_length|dna|Event_length|Center|Offset|Width|Residual|^BC[0-9]+$|^[A-Z]+[0-9]+DI$|bead|Dead|Live", colnames(data), ignore.case = T, invert = T)]
    data <- data[,which(toupper(colnames(data)) %in% toupper(input$selectmarkers))]
    data <- asinh(data / cofactor)
    if(length(unique(results$pheno_data$BATCH)) > 1){
      print("results$pheno_data batch:")
      print(head(results$pheno_data[match(files,results$pheno_data$SID),"BATCH"]))
      data <- apply(data, 2, function(x){x <- lm(x ~ results$pheno_data[match(files,results$pheno_data$SID),"BATCH"])$residual})
    }
    data <- data.frame(FILE = files, data)

    for (i in grep("FILE",colnames(data), ignore.case = T, invert = T)){
      data[,i] <- as.numeric(as.character(data[,i]))
    }
    
    counts <- data.frame(plyr::count(data$FILE))
    colnames(counts) <- c("SID","CELL_COUNT")
    results$pheno_data <- plyr::join(results$pheno_data, counts, by = "SID")
    
    sample_list <- unique(data$FILE)
    sample_list
    colnames(data)
    # marker_list <- as.character(colnames(data)[grep("FILE|population|tsne|SNEx|SNEy|PeakX|PeakY|time|length|DNA|^BC[0-9]+$|^[A-Z]+[0-9]+DI$|bead|L\\/D|Dead|Live", colnames(data), ignore.case = T, invert = T)])
    marker_list <- NULL
    marker_list <- input$selectmarkers
    results$selectedmarkers <- marker_list
    
    ggdf <- results$pheno_data
    ggdf <- ggdf[order(ggdf[,grep("Individual.*ID", colnames(ggdf), ignore.case = T)], decreasing = T),]
    ggdf[,grep("Sample.*ID", colnames(ggdf), ignore.case = T)] <- factor(ggdf[,grep("Sample.*ID", colnames(ggdf), ignore.case = T)], levels = unique(c(ggdf[,grep("Sample.*ID", colnames(ggdf), ignore.case = T)])))
    
    individual_colors <- gen_colors(color_conditions$general, length(unique(results$pheno_data$INDIVIDUAL_ID)))
    names(individual_colors) <- unique(results$pheno_data$INDIVIDUAL_ID)
    
    p1plots <- ggplot(ggdf, aes(x = SAMPLE_ID, y = CELL_COUNT, fill = INDIVIDUAL_ID, label = CELL_COUNT)) +
      geom_bar(stat="identity") +
      coord_flip() +
      theme_bw()+
      ylab("Cell Count")+
      scale_fill_manual(values = individual_colors) +
      guides(fill = guide_legend(reverse = TRUE))
    p1plots <- adjust_theme(p1plots)
    
    ggdf <- data.frame(SAMPLE_ID = results$pheno_data[match(data$FILE, results$pheno_data$SID),"SAMPLE_ID"],
                       data[,grep("FILE", colnames(data), ignore.case = T, invert = T)])
    
    ggdf <- reshape2::melt(ggdf, id.var = "SAMPLE_ID", value.name = "Expression", variable.name = "Marker")
    ggdf$GROUP <- results$pheno_data[match(ggdf$SAMPLE_ID, results$pheno_data$SAMPLE_ID),"GROUP"]
    
    group_colors <- gen_colors(color_conditions$bright[c(1,3:(length(color_conditions$bright)),2)], length(unique(results$pheno_data$GROUP)))
    names(group_colors) <- unique(results$pheno_data$GROUP)
    
    p2plots <- ggplot(ggdf, aes(x = Expression, y = SAMPLE_ID, color = GROUP, group = SAMPLE_ID)) +
      geom_density_ridges(alpha = 0.7) +
      facet_wrap(~ Marker, nrow = 4, scales = "free") +
      theme_bw() + ylab("Density") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 11),
            axis.text = element_text(size = 12),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            strip.text.x = element_text(size = 14),
            strip.background = element_blank(),
            legend.text = element_text(size=12),
            legend.title = element_text(size = 14),
            legend.key = element_rect(size = 6),
            legend.key.height = unit(1, "cm"),
            legend.key.width = unit(1, "cm")) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = group_colors)
    
    files <- unique(data$FILE)
    files
    plot_median <- NULL
    expr <- NULL
    
    for(i in 1:length(files)){
      current <- data[which(data$FILE == files[i]),]
      current <- current[,grep("FILE|peak|tsne|snex|sney|population",colnames(current), invert = T, ignore.case = T)]
      colnames(current)
      for (j in 1:ncol(current)){
        current[,j] <- as.numeric(as.character(current[,j]))
        expr <- c(expr,median(current[,j])) # current[,j] > 0
      }
      plot_median <- rbind(plot_median, expr)
      expr <- NULL
    }
    
    plot_median <- data.frame(t(plot_median))
    row.names(plot_median) <- marker_list
    colnames(plot_median) <- toupper(gsub("\\.fcs","",files,ignore.case = T))
    
    mds <- plotMDS(plot_median, plot = FALSE)
    pca_out <- prcomp(t(plot_median), center = TRUE, scale. = FALSE)
    ggdf <- data.frame(SID = colnames(plot_median), MDS1 = mds$x, MDS2 = mds$y, PC1 = pca_out$x[,1], PC2 = pca_out$x[,2])
    ggdf <- plyr::join(ggdf, results$pheno_data, by = "SID")
    
    p3plots <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = GROUP, label = SAMPLE_ID)) +
      geom_point(aes(size = CELL_COUNT), alpha = 0.8) +
      geom_label_repel(show.legend = F) +
      scale_size_continuous(range = c(6, 12))+
      # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
      theme_bw() +
      scale_color_manual(values = group_colors)
    p3plots <- adjust_theme(p3plots)
    
    p4plots <- ggplot(ggdf, aes(x = PC1, y = PC2, color = GROUP, label = SAMPLE_ID)) +
      geom_point(aes(size = CELL_COUNT), alpha = 0.8) +
      geom_label_repel(show.legend = F) +
      scale_size_continuous(range = c(6, 12))+
      # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
      theme_bw() +
      scale_color_manual(values = group_colors)
    p4plots <- adjust_theme(p4plots)
    
    print("plotmedian")
    print(head(plot_median))
    current <- plot_median
    current <- (scale((current)))
    current <- t(scale(t(current)))
    colnames(current) <- results$pheno_data[match(colnames(current), results$pheno_data$SID),"SAMPLE_ID"]
    p5plots <- complex_heatmap(as.matrix(current), col = jet2.col(n = 100, alpha = 1), legendtitle = "Z-Score")
    
    plotx <- melt(data.frame(Marker = row.names(plot_median), plot_median))
    colnames(plotx) <- c("MARKER","SID","ARCSINH_EXPRESSION")
    plotx <- cbind(plotx, results$pheno_data[match(gsub("\\.|\\s+","-",plotx$SID),gsub("\\.|\\s+","-",results$pheno_data$SID)),c("SAMPLE_ID","GROUP","BATCH")])
    
    p6plots <- ggboxplot(plotx, x = "GROUP", y = "ARCSINH_EXPRESSION", fill = "GROUP", palette = "jco", add = "jitter", repel = TRUE)
    p6plots <- p6plots + facet_wrap(~ MARKER, scales = "free", ncol = 6) + stat_compare_means() + theme_classic() +
      scale_fill_manual(values = group_colors)+
      scale_color_manual(values = results$sample_colors)
    p6plots <- adjust_theme(p6plots)
    # print(p6plots)
    
    # current <- compare_means(ARCSINH_EXPRESSION ~ GROUP, data = plotx, group.by = "MARKER")
    # current <- current[order(current$p.adj, decreasing = F),grep("\\.y\\.", colnames(current), ignore.case = T,invert = T)]
    p7data <- as.dendrogram(hclust(as.dist(1-cor((current)))))
    
    NRS <- function(x, ncomp = 3){
      pr <- prcomp(x, center = TRUE, scale. = FALSE)
      score <- rowSums(outer(rep(1, ncol(x)),pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
      return(score)
    }
    
    nrs_sample <- NULL
    current <- data[,grep("FILE|SNE|PEak|population",colnames(data), ignore.case = T, invert = T)]
    for(i in 1:length(unique(results$pheno_data$SID))){
      nrs_sample <- rbind(nrs_sample, NRS(current[which(data$FILE == unique(results$pheno_data$SID)[i]),]))
    }
    
    rownames(nrs_sample) <- unique(results$pheno_data$SID)
    nrs_sample <- data.frame(nrs_sample)
    nrs <- colMeans(nrs_sample, na.rm = TRUE)
    markers_ord <- names(sort(nrs, decreasing = TRUE))
    nrs_sample$SID <- rownames(nrs_sample)
    ggdf <- melt(nrs_sample, id.var = "SID",
                 value.name = "nrs", variable.name = "Markers")
    colnames(ggdf) <- c("SID","Markers","NRScore")
    ggdf$Markers <- factor(ggdf$Markers, levels = markers_ord)
    ggdf <- plyr::join(ggdf, results$pheno_data, by = "SID")
    
    p8plots <- ggplot(ggdf, aes(x = Markers, y = NRScore)) + 
      # geom_point(aes(color = SAMPLE_ID), alpha = 0.8,
      # position = position_jitter(width = 0.3, height = 0)) +
      geom_boxplot(aes(fill = GROUP), alpha = 0.8, outlier.color = NA) +
      theme_bw()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
      scale_color_manual(values = results$sample_colors)+
      scale_fill_manual(values = group_colors)
    p8plots <- adjust_theme(p8plots, xangle = 45, xsize = 12, hejust = 1, vejust = 1)
    
    files <- unique(data$FILE)
    data0 <- NULL
    
    for(i in 1:length(files)){
      data0[[i]] <- data[which(data$FILE == files[i]),grep("FILE|SNE|PEak|population",colnames(data), ignore.case = T, invert = T)]
    }
    
    for (i in 1:length(files)){
      meta <- data.frame(name=colnames(data0[[i]]),desc=colnames(data0[[i]]))
      meta$range <- apply(apply(data0[[i]],2,range),2,diff)
      meta$minRange <- apply(data0[[i]],2,min)
      meta$maxRange <- apply(data0[[i]],2,max)
      data0[[i]] <- new("flowFrame",exprs=as.matrix(data0[[i]]),parameters=AnnotatedDataFrame(meta))
    }
    
    removeModal()
    
    showModal(modalDialog("Estimating an optimal cluster..", footer=NULL))
    
    data_fs = as(data0,"flowSet")
    pData(data_fs)$name <- files
    som_input <- ReadInput(data_fs)
    set.seed(59)
    som <- BuildSOM(som_input)
    codes <- som$map$codes
    nmc <- 90
    mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 50,
                               pItem = 0.9, pFeature = 1, plot = NULL,
                               clusterAlg = "hc", innerLinkage = "complete", finalLinkage = "complete", distance = "euclidean", seed = 1234)
    Kvec = 2:nmc
    x1 = 0.1; x2 = 0.9
    PAC = rep(NA,length(Kvec)) 
    names(PAC) = paste("K=",Kvec,sep="")
    for(i in Kvec){
      M = mc[[i]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])
      PAC[i-1] = Fn(x2) - Fn(x1)
    }
    optK = Kvec[which.min(PAC)]
    PAC
    
    PAC <- data.frame(K = as.numeric(as.character(gsub("K=","",names(PAC)))), PAC = as.numeric(as.character(PAC)))
    ck1 <- elbow_point(PAC$K, PAC$PAC)
    ck1 <- ceiling(ck1$x)
    PAC$Difference <- NULL
    PAC[1,"Difference"] <- NA
    for (i in 2:nrow(PAC)){
      PAC[i,"Difference"] <- PAC[i-1,"PAC"] - PAC[i,"PAC"]
    }
    
    # cPAC <- PAC[(max(which(PAC$PAC > quantile(PAC$PAC, 0.95)))+1):nrow(PAC),]
    DeltaY <- diff(PAC$PAC)
    PAC_turn <- which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    ck2 <- min(PAC_turn[PAC_turn>5]) + 1
    if(ck1 != ck2){
      chosen_k <- min(ck1, ck2)
    }else{
      chosen_k <- ck1
    }
    
    PAC$Label <- ifelse(PAC$K == chosen_k, paste("Recommended K = ", chosen_k, sep = ""), "")
    
    p9plots <- ggplot(PAC, aes(x= K, y= PAC, label = Label)) + geom_line(colour = "grey") +
      ggtitle("Recommended K Number of Clusters Based on PAC Method")+
      theme_classic(base_size = 20)+ geom_text_repel(
        max.overlaps = Inf,force=1,
        point.padding = 0, # additional padding around each point
        min.segment.length = 0, # draw all line segments
        max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
        box.padding = 0.3, size = 12, colour = "red")+
      geom_point(colour = ifelse(PAC$K == chosen_k, "red", "grey"), size = ifelse(PAC$K == chosen_k, 5, 2))
    p9plots <- adjust_theme(p9plots)
    
    # results <- NULL
    results$data <- data
    results$p1plots <- p1plots
    results$p2plots <- p2plots
    results$p3plots <- p3plots
    results$p4plots <- p4plots
    results$p5plots <- p5plots
    results$p6plots <- p6plots
    results$p7data <- p7data
    results$p8plots <- p8plots
    results$p9plots <- p9plots
    results$individual_colors <- individual_colors
    results$group_colors <- group_colors
    results$recommend_k <- chosen_k
    results$mc <- mc
    results$som <- som
    removeModal()
    
  })
  
  observeEvent(input$submitclusters,{
    
    if(input$selectclusters == 0){
      showModal(modalDialog("Recommended clustering number will be used for the downstream analysis", footer=NULL, easyClose = T))
        updateTabsetPanel(session, "nav",selected = "section3")
        values$proceed3 <- 1
    }
    
    if(input$selectclusters != 0 & input$selectclusters < 2){
      showModal(modalDialog("Please either enter cluster 0 to use the recommended clustering number, or to enter a clustering number larger than 2", footer=NULL, easyClose = T))
      updateTabsetPanel(session, "nav",selected = "section3")
    }
    
    if(input$selectclusters > 1){
        updateTabsetPanel(session, "nav",selected = "section3")
        values$proceed3 <- 1
    }
  })
  
  observe({
    req(input$submitclusters)
    req(values$proceed3 == 1)
    showModal(modalDialog("Conducting analysis based on selected number of clusters..", footer=NULL))
    chosen_k <- NULL
    if(length(input$selectclusters) == 0){
      chosen_k <- results$recommend_k
    }else{
      chosen_k <- input$selectclusters
    }
    
    code_clustering1 <- results$mc[[chosen_k]]$consensusClass
    cell_clustering1 <- code_clustering1[results$som$map$mapping[,1]]
    results$data$population <- cell_clustering1

    plot_median <- NULL
    expr <- NULL
    cell_number <- NULL
    pop_list <- 1:chosen_k
    for(i in 1:length(pop_list)){
      current <- results$data[which(results$data$population == pop_list[i]),grep("FILE|SNE|PEak|population",colnames(results$data), ignore.case = T, invert = T)]
      cell_number <- c(cell_number,nrow(current))
      for (j in 1:ncol(current)){
        expr <- c(expr,median(current[,j])) # current[,j] > 0
      }
      plot_median <- rbind(plot_median, expr)
      expr <- NULL
    }
    
    row.names(plot_median) <- paste("Cluster_", pop_list, ":",cell_number,sep = "")
    colnames(plot_median) <- gsub("^.*?_","",results$selectedmarkers)
    plot_median <- data.frame(plot_median)
    dista <- hclust(as.dist(1-cor(t(plot_median))), method = "complete")
    plot_median[is.na(plot_median)] <- 0
    plot_median <- scale(plot_median)
    plot_median <- t(scale(t(plot_median)))
    plot_median <- as.matrix(plot_median)
    plot_median_csv <- plot_median[dista$order,]
    plot_median_csv <- plot_median_csv[,hclust(as.dist(1-cor((plot_median_csv))), method = "complete")$order]
    melt_plot_median <- melt(plot_median_csv)
    colnames(melt_plot_median) <- c("Cluster","Marker","Expression")
    melt_plot_median$Cluster_Order <- rep(1:nrow(plot_median_csv),ncol(plot_median_csv))
    melt_plot_median$Marker_Order <- rep(1:ncol(plot_median_csv),each = nrow(plot_median_csv))
    melt_plot_median$Size <- gsub(".*:(.*)","\\1",melt_plot_median$Cluster)
    melt_plot_median$Cluster <- gsub("Cluster_(.*):.*","\\1",melt_plot_median$Cluster)

    p10plots <- complex_heatmap(as.matrix(plot_median), col = jet2.col(n = 100, alpha = 1), legend_title = 'Z-Score')

    plot_median <- NULL
    expr <- NULL
    cell_number <- NULL
    pop_list <- 1:chosen_k
    files <- unique(results$data$FILE)
    for(i in 1:length(pop_list)){
      for(j in 1:length(files)){
        current <- results$data[which(results$data$population == pop_list[i] & results$data$FILE == files[j]),grep("FILE|SNE|PEak|population",colnames(results$data), ignore.case = T, invert = T)]
        if(nrow(current) > 0){
          cell_number <- c(cell_number,nrow(current))
          for (k in 1:ncol(current)){
            expr <- c(expr,median(current[,k])) # current[,k] > 0
          }
          plot_median <- rbind(plot_median, data.frame(FILE = files[j], CLUSTER = pop_list[i], MARKER = colnames(current), expr))
          expr <- NULL
        }
      }
    }
    
    colnames(plot_median) <- c("SID","CLUSTER","MARKER","ARCSINH_EXPRESSION")
    plot_median <- cbind(plot_median, results$pheno_data[match(gsub("\\.|\\s+","-",plot_median$SID),gsub("\\.|\\s+","-",results$pheno_data$SID)),c("SAMPLE_ID","GROUP","BATCH")])
    
    p11plots <- ggplot(plot_median) +
      geom_boxplot(aes(x = MARKER, y = ARCSINH_EXPRESSION, fill = GROUP), position = position_dodge(),
                   alpha = 1,
                   outlier.color = NA) +
      geom_point(aes(x = MARKER, y = ARCSINH_EXPRESSION), alpha = 0.8, size = 0.2) +
      facet_wrap(~ CLUSTER, scales = "free", ncol = 4) +
      theme_classic() +
      scale_fill_manual(values = results$group_colors)+
      scale_color_manual(values = results$sample_colors)
    p11plots <- adjust_theme(p11plots, xangle = 45, xsize = 8, title_size = 15, hejust = 1, vejust = 1)
    
    removeModal()
    
    showModal(modalDialog("Dimension reduction may take sometime depending on the size of the data and number of samples..", footer=NULL))
    
    print("Running PCA..")
    set.seed(59)
    pca_out <- prcomp(results$data[,grep("FILE|population|BARCODE", colnames(results$data), ignore.case = T, invert = T)], center = TRUE, scale. = FALSE)
    
    print("Running UMAP..")
    set.seed(59)
    umap_out <- umap::umap(results$data[,grep("FILE|population|BARCODE", colnames(results$data), ignore.case = T, invert = T)], pca = FALSE, perplexity = 30)
    
    plot_out <- data.frame(PC1 = pca_out$x[,1], PC2 = pca_out$x[,2],
                           UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2])
    
    plot_out$Cluster <- factor(results$data$population, levels = c(sort(as.numeric(as.character(unique(results$data$population))))))
    plot_out$FILE <- results$data$FILE
    plot_out$GROUP <- results$pheno_data[match(plot_out$FILE, results$pheno_data$SID),"GROUP"]
    plot_out$BATCH <- results$pheno_data[match(plot_out$FILE, results$pheno_data$SID),"BATCH"]
    
    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(results$data$population)))
    names(cluster_colors) <- c(sort(as.numeric(as.character(unique(results$data$population)))))
    
    p12plots <- ggplot(plot_out,  aes(x = UMAP1, y = UMAP2, color = GROUP)) +
      geom_point(size = 0.95, alpha = 1) +
      theme_bw() +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
      scale_color_manual(values = results$group_colors) +
      ggtitle(paste(results$project_name,": UMAP BY GROUP", sep = ""))
    p12plots <- adjust_theme(p12plots)
    
    p12plots2 <- ggplot(plot_out,  aes(x = UMAP1, y = UMAP2, color = BATCH)) +
      geom_point(size = 0.95, alpha = 1) +
      theme_bw() +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
      scale_color_manual(values = gen_colors(n = length(unique(plot_out$BATCH)))) +
      ggtitle(paste(results$project_name,": UMAP BY BATCH", sep = ""))
    p12plots2 <- adjust_theme(p12plots2)
    
    p12plots <- (p12plots | p12plots2)
    
    p13plots <- own_facet_scatter(plot_out, feature1 = "UMAP1", feature2 = "UMAP2", isfacet = T, legend_pos = "none",
                                  title = paste(results$project_name,": UMAP BY SAMPLE", sep = ""), col = results$sample_colors,ncol = 4,
                                  color_by = "FILE", group_by = "FILE", xlabel = "UMAP1", ylabel = "UMAP2", strip_size = 12)
    p13plots <- adjust_theme(p13plots, legend = "none", strip_size = 12)

    p14plots <- plot_bygroup(plot_out, x = "UMAP1", y = "UMAP2", group = "Cluster",
                             plot_title = paste(results$project_name,": UMAP BY CLUSTER", sep = ""),
                             col = cluster_colors, numeric = T, annot = T, point_size = 0.8, legendsize = 20, label_size = 8)
    p14plots <- adjust_theme(p14plots)
    
    p15plots <- plot_bygroup(plot_out, x = "PC1", y = "PC2", group = "Cluster",
                             plot_title = paste(results$project_name,": PCA BY CLUSTER", sep = ""),
                             col = cluster_colors, numeric = T, annot = F, point_size = 0.8, legendsize = 20, label_size = 8)
    p15plots <- adjust_theme(p15plots)
    
    marker_list <- unique(colnames(results$data[,grep("FILE|population", colnames(results$data), ignore.case = T, invert = T)]))
    melt_data <- melt(data.frame(results$data, UMAP1 = plot_out$UMAP1, UMAP2 = plot_out$UMAP2), 
                      measure.vars = colnames(results$data)[grep("FILE|population", colnames(results$data), ignore.case = T, invert = T)], id.vars = c("UMAP1","UMAP2", "FILE"))
    colnames(melt_data) <- c("UMAP1","UMAP2","FILE","Marker","Asinh_Expression")
    melt_data$GROUP <- results$pheno_data[match(melt_data$FILE, results$pheno_data$SID),"GROUP"]
    p16data <- melt_data
    
    sample_list <- sort(as.character(unique(results$data$FILE)))
    sample_list
    summary <- data.frame(population = sort(unique(results$data$population)))
    
    for (i in 1:length(sample_list)){
      current <- results$data[which(results$data$FILE == sample_list[i]),]
      summary <- plyr::join(summary,plyr::count(current, c("FILE","population"))[c("population","freq")],by = "population", type = "left")
      colnames(summary)[which(colnames(summary) == "freq")] <- unique(current$FILE)
    }
    summary[is.na(summary)] <- 0
    
    colnames(summary)[which(colnames(summary) == "population")] <- "Cluster"
    summary_percent = data.frame(summary)
    summary_percent$Total_Cells_Number <- rowSums(summary_percent[,which(!colnames(summary) == "Cluster")])
    
    disease_types <- as.character(unique(results$pheno_data$GROUP))
    
    for (i in 1:length(disease_types)){
      mm <- c("Cluster",as.character(results$pheno_data[match(colnames(summary)[grep("Cluster",colnames(summary), invert = T, ignore.case = T)],as.character(results$pheno_data$SID)),"GROUP"]))
      if(length(grep(as.character(disease_types[i]),mm, ignore.case = T)) > 1){
        current_rs <- rowSums(summary[,grep(as.character(disease_types[i]),mm, ignore.case = T)])
      } else{
        current_rs <- as.numeric(as.character(summary[,grep(as.character(disease_types[i]),mm, ignore.case = T)]))
      }
      summary_percent[,ncol(summary_percent)+1] <- current_rs
      colnames(summary_percent)[ncol(summary_percent)] <- paste("GROUP_",disease_types[i],"_Total_Cells", sep = "")
      summary_percent[,ncol(summary_percent)+1] <- (current_rs/summary_percent[,"Total_Cells_Number"])*100
      colnames(summary_percent)[ncol(summary_percent)] <- paste("GROUP_",disease_types[i],"_Percent", sep = "")
    }
    
    for(l in grep("Cluster|Total_Cells|Strata_|Percent",colnames(summary_percent), ignore.case = TRUE, invert = TRUE)){
      summary_percent[,l] <- summary_percent[,l]/summary_percent[,"Total_Cells_Number"]
    }
    
    melt_summary <- melt(summary_percent[,grep("Percent|Total_Cells|Strata_",colnames(summary_percent), ignore.case = TRUE, invert = T)], id.vars = c("Cluster"))
    colnames(melt_summary) <- c("Cluster","SID","Proportion")
    melt_summary <- cbind(melt_summary, results$pheno_data[match(gsub("\\.|\\s+","\\-",melt_summary$SID), gsub("\\.|\\s+","\\-",results$pheno_data$SID)),grep("SID", colnames(results$pheno_data), ignore.case = T, invert = T)])
    melt_summary$Cluster <- factor(melt_summary$Cluster, levels = sort(unique(as.numeric(as.character(melt_summary$Cluster)))))

    p17plots <- ggplot(melt_summary, aes(Cluster, Proportion, fill = SAMPLE_ID))+
      geom_bar(stat="identity", alpha=0.8)+
      coord_polar()+
      scale_fill_viridis(discrete=TRUE)+
      ggtitle(paste(results$project_name,"\nFrequency of Samples in Each Cluster", sep = ""))+
      theme_bw(base_size = 28)+
      theme(plot.margin = unit(c(3,3,3,3), "cm"),
            plot.title = element_text(size=25, face = "bold", hjust = 0.5),
            legend.title=element_text(size=20, face = "bold"), 
            legend.text=element_text(size=18),
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
    
    p17bplots <- ggplot(melt_summary,aes(Cluster,Proportion,fill=SAMPLE_ID))+ 
      geom_bar(stat = "identity", position = "fill", color = "black", size = 1)+ #, color = "black", size = 1
      theme_classic()+
      scale_fill_viridis(discrete=TRUE)+
      guides(color=guide_legend(title="SAMPLES", ncol = 2),fill = guide_legend(title="SAMPLES", ncol = 2))
    p17bplots <- adjust_theme(p17bplots, xangle = 45, hejust = 1, vejust = 1, strip_size = 1)
    
    p18plots <- ggplot(melt_summary,aes(SAMPLE_ID,Proportion,fill=Cluster))+ 
      geom_bar(stat = "identity", position = "fill", color = "black", size = 1)+ #, color = "black", size = 1
      theme_classic()+
      scale_fill_manual(values = cluster_colors)+
      guides(color=guide_legend(title="CLUSTERS", ncol = 2),fill = guide_legend(title="CLUSTERS", ncol = 2))
    p18plots <- adjust_theme(p18plots, xangle = 45, hejust = 1, vejust = 1, strip_size = 1)
    
    codes <- NULL
    codes <- results$som$map$codes
    code_sizes <- table(factor(results$som$map$mapping[,1], levels = 1:nrow(codes)))
    code_sizes <- as.numeric(as.character(code_sizes))
    
    removeModal()
    
    showModal(modalDialog("Running dimension reduction on meta clustering..", footer=NULL))
    print("Running UMAP..")
    set.seed(59)
    umap_out <- umap(codes)
    # tsne_out <- Rtsne(codes)
    pca_out <- prcomp(codes, center = TRUE, scale. = FALSE)
    codes_dr <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], 
                           # tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                           PC1 = pca_out$x[, 1], PC2 = pca_out$x[, 2])
    codes_dr$Cluster <- factor(code_clustering1)
    codes_dr$Size <- code_sizes
    codes_dr$SOM_ID <- as.numeric(as.character(row.names(codes_dr)))
    
    # som_data <- data
    results$data$som_population <- results$som$map$mapping[,1]
    som_result <- data.frame(SOM_ID = sort(unique(results$data$som_population)))
    som_result$SOM_Phenotype <- "UNDEFINED"
    
    # choose top n to display
    top_n_markers <- 5
    # markers_to_note <- c("CD4","CD8","CD56","CD45RA","CD45RO")
    
    current0 = apply(results$data[,which(toupper(colnames(results$data)) %in% toupper(marker_list))], 2, list)
    current0 = lapply(current0, function(x){x = (as.numeric(unlist(x))[as.numeric(unlist(x))>0])})
    median_c0 <- as.numeric(lapply(current0, median))
    
    som_result$Cluster <- codes_dr[match(som_result$SOM_ID, codes_dr$SOM_ID),"code_clustering1"]
    
    node_result <- data.frame(Cluster_ID = sort(unique(results$data$population)))
    node_result$Cluster_Phenotype <- "UNDEFINED"
    
    for(i in 1:nrow(node_result)){
      node_name3 <- NULL
      current_node <- as.numeric(as.character(node_result[i,"Cluster_ID"]))
      current3 = apply(results$data[results$data$population==current_node,grep(paste(marker_list, collapse = "|"), colnames(results$data), ignore.case = TRUE)], 2, list)
      current3 = lapply(current3, function(x){x = (as.numeric(unlist(x))[as.numeric(unlist(x))>1])})
      length_c3 = as.numeric(unlist(lapply(current3, length)))
      median_c3 <- as.numeric(lapply(current3, median))
      median_c3c0 <- (median_c3+1)/(median_c0+1)
      node_list3 = marker_list[ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE)][order(median_c3c0[which(ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE))], decreasing = T)]
      if (length(node_list3) == 0){
        node_name3 = "UNDEFINED"
      }else{
        if(length(node_list3) > top_n_markers){node_list3 <- node_list3[1:top_n_markers]}
        node_name3 <- paste(node_name3,paste(node_list3, collapse = "+"), "+", sep = "")
      }
      node_result[i,"Cluster_Phenotype"] <- node_name3
      node_list3 <- NULL
    }
    
    codes_dr$Cluster_Phenotype <- node_result[match(codes_dr$Cluster, node_result$Cluster_ID),"Cluster_Phenotype"]
    colnames(codes_dr)[which(colnames(codes_dr) == "Cluster")] <- "Cluster"
    
    for(i in 1:nrow(som_result)){
      node_name3 <- NULL
      current_node <- as.numeric(as.character(som_result[i,"SOM_ID"]))
      current3 = apply(results$data[results$data$som_population==current_node,grep(paste(marker_list, collapse = "|"), colnames(results$data), ignore.case = TRUE)], 2, list)
      current3 = lapply(current3, function(x){x = (as.numeric(unlist(x))[as.numeric(unlist(x))>0])})
      length_c3 = as.numeric(unlist(lapply(current3, length)))
      median_c3 <- as.numeric(lapply(current3, median))
      median_c3c0 <- (median_c3+1)/(median_c0+1)
      node_list3 = marker_list[ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE)][order(median_c3c0[which(ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE))], decreasing = T)]
      if (length(node_list3) == 0){
        node_name3 = "UNDEFINED"
      }else{
        if(length(node_list3) > top_n_markers){node_list3 <- node_list3[1:top_n_markers]}
        node_name3 <- paste(node_name3,paste(node_list3, collapse = "+"), "+", sep = "")
        
      }
      som_result[i,"SOM_Phenotype"] <- node_name3
      node_list3 <- NULL
      
    }
    
    codes_dr$SOM_Phenotype <- som_result[match(codes_dr$SOM_ID, som_result$SOM_ID),"SOM_Phenotype"]
    # For user to download their SOM Phenotypes
    # codes_dr_simp <- codes_dr[,c("SOM_ID","SOM_Phenotype","Cluster","Cluster_Phenotype","Size")]
    # colnames(codes_dr_simp) <- c("SOM_ID","SOM_Phenotype","NODE_ID","NODE_Phenotype","Cluster_Size")
    codes_dr$SOM_ID <- factor(codes_dr$SOM_ID, levels = sort(unique(as.numeric(as.character(codes_dr$SOM_ID)))))
    colnames(codes_dr)[grep("code_clustering", colnames(codes_dr), ignore.case = T)] <- "Cluster"
    codes_dr$Cluster <- factor(codes_dr$Cluster, levels = sort(unique(as.numeric(as.character(codes_dr$Cluster)))))
    colnames(codes_dr)[which(colnames(codes_dr) == "SOM_ID")] <- "Meta_Cluster"
    colnames(codes_dr)[which(colnames(codes_dr) == "Cluster")] <- "Final_Cluster"
    colnames(codes_dr)[which(colnames(codes_dr) == "SOM_Phenotype")] <- "Meta_Cluster_Phenotype"
    colnames(codes_dr)[which(colnames(codes_dr) == "Cluster_Phenotype")] <- "Final_Cluster_Phenotype"
    codes_dr$Label <- paste("Meta_Cluster:",codes_dr$Meta_Cluster, sep = ":")
    p20data <- unique(codes_dr[,grep("UMAP|tSNE|PC[0-9]+|PCA|Label", colnames(codes_dr), ignore.case = T, invert = T)])
    p20data <- p20data[,c("Meta_Cluster","Size","Meta_Cluster_Phenotype","Final_Cluster","Final_Cluster_Phenotype")]
    # p19plots <- ggplot(codes_dr, aes(x = tSNE1, y = tSNE2,color = Cluster, size = Size, label = Label)) +
    #   geom_point(alpha = 0.7) + theme_classic() +
    #   scale_point_size_continuous(range = c(9,18))+
    #   geom_label_repel(box.padding = 0.2, max.overlaps = Inf, size = 2.2, show.legend = F)+
    #   scale_color_manual(values = cluster_colors) +
    #   theme(legend.position = "right",
    #         plot.title = element_text(size=18, face = "bold", hjust = 0.5))+
    #   ggtitle(paste(project_name,"\nSOM Nodes tSNE", sep = ""))
    # p19plots <- adjust_theme(p19plots)
    # print("codes_dr$Label:")
    # print(table(codes_dr$Label))
    
    p21plots <- ggplot(codes_dr, aes(x = UMAP1, y = UMAP2,color = Final_Cluster, size = Size, label = Label)) +
      geom_point(alpha = 0.7) + theme_classic() +
      scale_point_size_continuous(range = c(9,18))+
      geom_text_repel(box.padding = 0.2, max.overlaps = Inf, show.legend = F)+
      scale_color_manual(values = cluster_colors) +
      theme(legend.position = "right",
            plot.title = element_text(size=18, face = "bold", hjust = 0.5))+
      ggtitle(paste(results$project_name,"\nUMAP OF META-CLUSTERS COLORED BY FINAL CLUSTERS", sep = ""))
    p21plots <- adjust_theme(p21plots)
    
    p22plots <- ggplot(codes_dr, aes(x = PC1, y = PC2, color = Final_Cluster, size = Size, label = Label)) +
      geom_point(alpha = 0.7) + theme_classic() +
      scale_point_size_continuous(range = c(9,18))+
      geom_text_repel(box.padding = 0.2, max.overlaps = Inf, show.legend = F)+
      scale_color_manual(values = cluster_colors) +
      theme(legend.position = "right",
            plot.title = element_text(size=18, face = "bold", hjust = 0.5))+
      ggtitle(paste(results$project_name,"\nPCA OF META-CLUSTERS COLORED BY FINAL CLUSTERS", sep = ""))
    p22plots <- adjust_theme(p22plots)

    # results <- NULL
    # results$data <- data
    # results$p1plots <- p1plots
    # results$p2plots <- p2plots
    # results$p3plots <- p3plots
    # results$p4plots <- p4plots
    # results$p5plots <- p5plots
    # results$p6plots <- p6plots
    # results$p7data <- p7data
    # results$p8plots <- p8plots
    # results$p9plots <- p9plots
    results$p10plots <- p10plots
    results$p11plots <- p11plots
    results$p11data <- plot_median
    results$p12plots <- p12plots
    results$p13plots <- p13plots
    results$p14plots <- p14plots
    results$p15plots <- p15plots
    results$p16data <- p16data
    results$p17plots <- p17plots
    results$p17bplots <- p17bplots
    results$p18plots <- p18plots
    # results$p19plots <- p19plots
    results$p20data <- p20data
    results$p21plots <- p21plots
    results$p22plots <- p22plots
    results$cluster_colors <- cluster_colors
    #######################################################################################################################################
    # saveRDS(results,"results.RDS")
    groups <- sort(unique(results$p16data$GROUP))
    clusters <- sort(unique(results$p20data$Final_Cluster))
    # annot_names <- sort(unique(results$data$FILE))
    removeModal()
    
    updateSelectInput(session, inputId = 'p16id', label = 'Choose a group to display', choices = groups, selected = groups[1])
    updateSelectInput(session, inputId = 'p23id', label = 'Choose a cluster to display', choices = clusters, selected = clusters[1])
    updateSelectInput(session, inputId = 'p11bid', label = 'Choose a cluster to display', choices = clusters, selected = clusters[1])
    
    print("Completed!")
    return(results)
  })
  
  output$p1plot <- renderPlot({ #renderPlotly
    showModal(modalDialog("Plotting figure 1..", footer=NULL))
    p <- results[['p1plots']]
    return(p)
    removeModal()
  }, height = 1000, width = 1400)
  
  observeEvent(input$p1plot, {
    screenshot(id="p1plot", filename = paste("1SCA_CyTOF_SAMPLE_CELL_SUMMARY_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p2plot <- renderPlot({
      showModal(modalDialog("Plotting figure 2..", footer=NULL))
      p <- results[['p2plots']]
      return(p)
      removeModal()
  }, height = 1000, width = 1400)
  
  observeEvent(input$p2plot, {
    screenshot(id="p2plot", filename = paste("2SCA_CyTOF_MARKER_DENSITIES_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p3plot <- renderPlot({
      showModal(modalDialog("Plotting figure 3..", footer=NULL))
      p <- results[['p3plots']]
      removeModal()
      return(p)
  }, height = 800, width = 1400)
  
  observeEvent(input$p3plot, {
    screenshot(id="p3plot", filename = paste("3SCA_CyTOF_MDS_SAMPLE_DISTANCE_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p4plot <- renderPlot({
      showModal(modalDialog("Plotting figure 4..", footer=NULL))
      p <- results[['p4plots']]
      removeModal()
      return(p)
  }, height = 800, width = 1400)
  
  observeEvent(input$p4plot, {
    screenshot(id="p4plot", filename = paste("4SCA_CyTOF_PCA_SAMPLE_DISTANCE_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p5plot <- renderPlot({
    showModal(modalDialog("Plotting figure 5..", footer=NULL))
    p <- results[['p5plots']]
    removeModal()
    return(p)
  }, height = 800, width = 1400)
  
  observeEvent(input$p5plot, {
    screenshot(id="p5plot", filename = paste("5SCA_CyTOF_HEATMAP_SAMPLE_DISTANCE_ARCSINH_MEDIAN_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p6plot <- renderPlot({
    showModal(modalDialog("Plotting figure 6..", footer=NULL))
    p <- results[['p6plots']]
    removeModal()
    return(p)
  }, height = 1400, width = 1400)
  
  observeEvent(input$p6plot, {
    screenshot(id="p6plot", filename = paste("6SCA_CyTOF_BOXPLOT_MARKERS_BY_GROUP_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p7plot <- renderPlot({
    showModal(modalDialog("Plotting figure 7..", footer=NULL))
    par(mar=c(3,4,1,6))
    print(plot(results[['p7data']], horiz = TRUE))
    removeModal()
    return(p)
  }, height = 1000, width = 1000)
  
  observeEvent(input$p7plot, {
    screenshot(id="p7plot", filename = paste("7SCA_CyTOF_DENDROGRAM_SAMPLES_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p8plot <- renderPlot({
    showModal(modalDialog("Plotting figure 8..", footer=NULL))
    p <- results[['p8plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p8plot, {
    screenshot(id="p8plot", filename = paste("8SCA_CyTOF_NRS_RANKING_ALL_MARKERS_",results[['project_name']], sep = ""), scale = 2)
  })

output$p9plot <- renderPlot({
  showModal(modalDialog("Plotting figure 9..", footer=NULL))
  p <- results[['p9plots']]
  removeModal()
  return(p)
}, height = 700, width = 1400)

observeEvent(input$p9plot, {
  screenshot(id="p9plot", filename = paste("9SCA_CyTOF_PAC_ELBOW_PLOT_CHOSEN_CLUSTER_NUMBER_",results[['project_name']], sep = ""), scale = 2)
})

output$p10plot <- renderPlot({
  showModal(modalDialog("Plotting figure 10..", footer=NULL))
  p <- results[['p10plots']]
  removeModal()
  return(p)
}, height = 700, width = 1400)

observeEvent(input$p10plot, {
  screenshot(id="p10plot", filename = paste("10SCA_CyTOF_HEATMAP_CLUSTERS_ARCSINH_MEDIAN_",results[['project_name']], sep = ""), scale = 2)
})

output$p11plot <- renderPlot({
  showModal(modalDialog("Plotting figure 11..", footer=NULL))
  p <- results[['p11plots']]
  removeModal()
  return(p)
}, height = 1400, width = 1400)

observeEvent(input$p11plot, {
  screenshot(id="p11plot", filename = paste("11SCA_CyTOF_BOXPLOT_CLUSTERS_BY_GROUP_",results[['project_name']], sep = ""), scale = 2)
})

output$p11bplot <- renderPlot({
  showModal(modalDialog("Plotting figure 11 clusterwise plots..", footer=NULL))
  plotx <- results[['p11data']]
  plotx <- plotx[which(plotx$CLUSTER == input$p11bid),]
  p <- ggboxplot(plotx, x = "GROUP", y = "ARCSINH_EXPRESSION", fill = "GROUP", palette = "jco", add = "jitter", repel = TRUE)
  p <- p + facet_wrap(~ MARKER, scales = "free", ncol = 6) + stat_compare_means() + theme_classic() +
    scale_fill_manual(values = results[['group_colors']])+
    scale_color_manual(values = results[['sample_colors']])
  p <- adjust_theme(p, xsize = 12)
  removeModal()
  return(p)
}, height = 900, width = 1200)

observeEvent(input$p11plot, {
  screenshot(id="p11plot", filename = paste("11SCA_CyTOF_BOXPLOT_CLUSTERS_BY_GROUP_",results[['project_name']], sep = ""), scale = 2)
})

output$p12plot <- renderPlot({
  showModal(modalDialog("Plotting figure 12..", footer=NULL))
  p <- results[['p12plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p12plot, {
  screenshot(id="p12plot", filename = paste("12SCA_CyTOF_UMAP_BY_GROUP_",results[['project_name']], sep = ""), scale = 2)
})

output$p13plot <- renderPlot({
  showModal(modalDialog("Plotting figure 13..", footer=NULL))
  p <- results[['p13plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p13plot, {
  screenshot(id="p13plot", filename = paste("13SCA_CyTOF_UMAP_BY_SAMPLE_",results[['project_name']], sep = ""), scale = 2)
})

output$p14plot <- renderPlot({
  showModal(modalDialog("Plotting figure 14..", footer=NULL))
  p <- results[['p14plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p14plot, {
  screenshot(id="p14plot", filename = paste("14SCA_CyTOF_UMAP_BY_CLUSTER_",results[['project_name']], sep = ""), scale = 2)
})

output$p15plot <- renderPlot({
  showModal(modalDialog("Plotting figure 15..", footer=NULL))
  p <- results[['p15plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p15plot, {
  screenshot(id="p15plot", filename = paste("15SCA_CyTOF_PCA_BY_CLUSTER_",results[['project_name']], sep = ""), scale = 2)
})

output$p16plot <- renderPlot({
if(input$p16id != ""){
  print("input p16id:")
  print(input$p16id)
  showModal(modalDialog("Plotting figure 16..", footer=NULL))
  plotx <- results[['p16data']][which(results[['p16data']]$GROUP == input$p16id),]
  p <- ggplot(plotx,  aes(x = UMAP1, y = UMAP2, color = Asinh_Expression)) +
    facet_wrap(~Marker, ncol = 6) +
    geom_point(size = 0.8) + theme_classic() + ggtitle(input$p16id) +
    scale_color_gradientn(colours = colorRampPalette(color_conditions$BlueYellowRed)(100))
  p <- adjust_theme(p)
  removeModal()
  return(p)
}
}, height = 1400, width = 1000)

observeEvent(input$p16plot, {
  screenshot(id="p16plot", filename = paste("16SCA_CyTOF_UMAP_ARCSINH_EXPRESSION_MARKER_",input$p16id, sep = ""), scale = 2)
})

output$p17plot <- renderPlot({
  showModal(modalDialog("Plotting figure 17..", footer=NULL))
  p <- results[['p17plots']]
  removeModal()
  return(p)
}, height = 1400, width = 1400)

observeEvent(input$p17plot, {
  screenshot(id="p17plot", filename = paste("17SCA_CyTOF_PROPORTION_OF_SAMPLES_IN_EACH_CLUSTER_PIE_",results[['project_name']], sep = ""), scale = 2)
})

output$p17bplot <- renderPlot({
  showModal(modalDialog("Plotting figure 17 stacked bar plot..", footer=NULL))
  p <- results[['p17bplots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p17bplot, {
  screenshot(id="p17bplot", filename = paste("17SCA_CyTOF_PROPORTION_OF_SAMPLES_IN_EACH_CLUSTER_BARPLOT_",results[['project_name']], sep = ""), scale = 2)
})

output$p18plot <- renderPlot({
  showModal(modalDialog("Plotting figure 18..", footer=NULL))
  p <- results[['p18plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p18plot, {
  screenshot(id="p18plot", filename = paste("18SCA_CyTOF_PROPORTION_OF_CLUSTERS_IN_EACH_SAMPLE_",results[['project_name']], sep = ""), scale = 2)
})

# output$p19plot <- renderPlot({
#   showModal(modalDialog("Plotting figure 19..", footer=NULL))
#   p <- results[['p19plots']]
#   removeModal()
#   return(p)
# }, height = 800, width = 1400)
# 
# observeEvent(input$p19plot, {
#   screenshot(id="p19plot", filename = paste("19SCA_CyTOF_tSNE_SOM_NODES_",results[['project_name']], sep = ""), scale = 2)
# })

output$p20table <- renderDT(results[['p20data']],
                            filter = "top",
                            style="bootstrap",
                            rownames = F,
                            options = list(pageLength = 10))

output$p21plot <- renderPlot({
  showModal(modalDialog("Plotting figure 21..", footer=NULL))
  p <- results[['p21plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p21plot, {
  screenshot(id="p21plot", filename = paste("21SCA_CyTOF_UMAP_SOM_NODES_",results[['project_name']], sep = ""), scale = 2)
})

output$p22plot <- renderPlot({
  showModal(modalDialog("Plotting figure 22..", footer=NULL))
  p <- results[['p22plots']]
  removeModal()
  return(p)
}, height = 800, width = 1400)

observeEvent(input$p22plot, {
  screenshot(id="p22plot", filename = paste("22SCA_CyTOF_PCA_SOM_NODES_",results[['project_name']], sep = ""), scale = 2)
})

output$p23plot <- renderPlot({
  if(input$p23id != ""){
    print("input p23id:")
    print(input$p23id)
    showModal(modalDialog("Plotting figure 23..", footer=NULL))
    data <- results[['data']]
    data_1 <- data[which(data$population == input$p23id),]
    nodes <- sort(unique(data$population))
    par(mfrow=c(ceiling(length(nodes)/4),4))
    par(mar=c(4,2,2,1))
    dis = rep(0,length(colnames(data)[grep("FILE|SAMPLE.*ID|population|tSNE",colnames(data),ignore.case = T,invert = T)]))
    
    for(j in grep("FILE|SAMPLE.*ID|population|tSNE",colnames(data),ignore.case = T,invert = T)){
      x=data[data[,j]>0,j]
      y=data_1[data_1[,j]>0,j]
      if (length(x)>3)
      {
        dis[j-1]=median(y)/median(x)
      }
    }
    
    for(j in grep("FILE|SAMPLE.*ID|population|tSNE",colnames(data),ignore.case = T,invert = T))
    {
      pos=order(-dis)[j-1]+1
      x=data[data[,pos]>0,pos]
      y=data_1[data_1[,pos]>0,pos]
      
      if (length(x)>3)
      {
        plot(density(x),xlab='',ylab='Density',main=colnames(data)[pos],
             xlim=c(-1,10),ylim=c(0,3))
        abline(v=median(x),col=1)
        if(length(y)>3) lines(density(y),col='red')
        abline(v=median(y),col=2)
        text(7,2.9,paste('All(',nrow(data),')=',
                         as.character(round(median(x), digits = 1)),'',sep=''),col=1)
        text(7,2.5,paste('Cluster ',input$p23id,'(',dim(data_1)[1],')=',
                         as.character(round(median(y), digits = 1)),'',sep=''),col=2)
      }
    }
    removeModal()
  }
}, width = 900, height = 1200)

observeEvent(input$p23plot, {
  screenshot(id="p23plot", filename = paste("23SCA_CyTOF_MARKER_DENSITY_PLOTS_CLUSTER_",input$p23id, sep = ""), scale = 2)
})

})


