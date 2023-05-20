###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: scCNV Analysis Pipeline
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
library("GenomicRanges")
library("pheatmap")
library("cowplot")
library("ComplexHeatmap")
library("patchwork")
library("jamba")
library("splicejam")
library("ape")
library("data.table")
library("intervalaverage")
library("ggplot2")
library("RColorBrewer")
library("ggridges")
library("adegenet")
library("reshape2")
library("ggtree")
library("data.table")
library("Seurat")
library("fastcluster")

source("DB/SCA_scCNV_RShiny_Functions_V1.0.0.R")
color_conditions <- color_ini()

ctime <- format(Sys.time(), format = "%Y%m%d%H%M%S", tz = "Europe/Stockholm")

ploidy_levels <- c("0","1","2","3","4","5","6","7","8","9",">=10")

color_heatmap <- c("#306192","#548292","white","#7EB0C0","#9BBAA0","#B5C07F","#C8C667","#DEBD42","#D8AC3A","#D96829","#DE3C22")
names(color_heatmap) <- ploidy_levels

cell_ploidy <- c("diploid", "non-diploid", "noisy")
cell_ploidy_colors <- c("#5F75A5", "#B0A062", "#D3B6A5")
names(cell_ploidy_colors) <- cell_ploidy

size_threshold <- 2000000
confidence_threshold <- 20

example1 <- "DB/SCA_scCNV_Example_From_10X.zip"
example2 <- "DB/SCA_scCNV_Metadata_Example.csv"

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
    filename = "SCA_scCNV_Example_From_10X.zip",
    content = function(file) {
      file.copy(example1, file)
    },
    contentType = "application/zip"
  )
  
  output$downloadExample2 <- downloadHandler(
    filename = "SCA_scCNV_Metadata_Example.csv",
    content = function(file) {
      file.copy(example2, file)
    },
    contentType = "text/csv"
  )
  
  data <- reactive({
    req(input$files)
    showModal(modalDialog("Processing..depending on the number of samples as well as file sizes, time of processing may varies. Please wait patiently for results to be delivered.", footer=NULL))
    inFile <- input$files
    pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "CNV", isDir = T)
    print(input$files)
    print(input$files$datapath)
    project_name <- ifelse(input$project_name == "","My_Project",gsub("\\s+|\\(|\\)|-|\\/|\\?","",input$project_name))
    print(project_name)
    
    cdir <- paste(getwd(),"/",project_name,ctime, "/",sep = "")
    print("cdir:")
    print(cdir)
    dir.create(cdir)
    print("inFile$datapath:")
    print(inFile$datapath)
    print("Unzipping..")
    sample_files <- unzip(inFile$datapath, exdir = cdir)
    sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
    print("sample files:")
    print(sample_files)
    input_dir <- gsub("(.*\\/).*$","\\1",sample_files[1], ignore.case = T)
    print("input_dir:")
    print(input_dir)
    
    ##################################################################################################
    data <- NULL
    data_current <- NULL
    data_summary <- NULL
    data_cell_stats <- NULL
    data_mappable <- NULL
    data_group <- NULL
    hc_summary <- NULL
    
    removeModal()
    
    for(i in 1:nrow(pheno_data)){
      print(pheno_data[i,"SAMPLE_ID"])
      showModal(modalDialog(paste("Processing sample: ", pheno_data[i,"SAMPLE_ID"], "..",sep = ""), footer=NULL))
      current_dir <- list.files(path = input_dir, pattern = pheno_data[i,grep("FILE", colnames(pheno_data), ignore.case = T)], ignore.case = T, full.names = T)
      if(length(current_dir) > 0){
        current_bed <- current_dir[grep(current_dir, pattern = "\\cnv.*calls.*bed$", ignore.case = T)][1]
        print("current_bed")
        print(current_bed)
        print("reading in bed..")
        current <- read.table(current_bed)
        print("finished reading in bed..")
        colnames(current) <- c("chrom","start","end","id","copy_number","event_confidence")
        current_summary <- current_dir[grep(current_dir, pattern = "summary_metrics.*\\.csv$", ignore.case = T)][1]
        current_summary <- read.csv(current_summary)
        current_mappable <- current_dir[grep(current_dir, pattern = "\\mappable_regions.*bed$", ignore.case = T)][1]
        print("current_mappable:")
        print(current_mappable)
        print("reading in mappable..")
        current_mappable <- read.table(current_mappable)
        colnames(current_mappable) <- c("chrom","start","end")
        data_mappable <- rbind(data_mappable, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current_mappable))
        current_group <- current[which(!((current$id %in% current_summary$cell_id) & (current$id <= max(current_summary$cell_id)))),]
        data_group <- rbind(data_group, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current_group))
        current$Size <- as.numeric(as.character(current$end)) - as.numeric(as.character(current$start))
        # current <- current[which(current$Size >= size_threshold & current$event_confidence > confidence_threshold),]
        current <- split(current, current$chrom)
        current <- lapply(current, function(x){x <- x[order(x$start, decreasing = F),]})
        current <- do.call(rbind.data.frame,current)
        current <- current[which(current$id %in% current_summary$cell_id),]
        current$Ploidy <- current$copy_number
        current[which(current$Ploidy >=10),"Ploidy"] <- ">=10"
        current$Ploidy <- factor(current$Ploidy, levels = ploidy_levels)
        current$Region_type <- ifelse(current$copy_number == 2, "diploid","non-diploid")
        current[current_summary[match(current$id,current_summary$cell_id),"is_noisy"] == 1,"Region_type"] <- "noisy"
        chr_num <- sort(as.numeric(as.character(unique(current$chrom)[grep("X|Y",unique(current$chrom), ignore.case = T, invert = T)])))
        chr_num <- c(as.character(chr_num), unique(current$chrom)[grep("X|Y",unique(current$chrom), ignore.case = T)])
        
        temp <- split(current, current$chrom)
        temp <- lapply(temp, function(x){
          x <- data.table(x)
          colnames(x)[grep("start", colnames(x), ignore.case = T)] <- "start_old"
          colnames(x)[grep("end", colnames(x), ignore.case = T)] <- "end_old"
          x <- isolateoverlaps(x,interval_vars = c("start_old","end_old"),group_vars = "id")
        })
        
        paste("Finnished isolating overlaps ..")
        current <- do.call(rbind.data.frame,temp)
        current$range <- paste(current$chrom,":",current$start,"-",current$end, sep = "")
        duplicated_data <- NULL
        duplicated_data <- current[!(!(duplicated(current[,c("id","range")]) | duplicated(current[,c("id","range")], fromLast = TRUE))),]
        if(nrow(duplicated_data) > 0){
          duplicated_data <- split(duplicated_data, duplicated_data$id)
          duplicated_data <- lapply(duplicated_data, function(x){x <- split(x, x$range)})
          duplicated_data <- lapply(duplicated_data, function(x){lapply(x, function(y){
            y <- y[which.max(y$event_confidence),]
          })})
          duplicated_data <- lapply(duplicated_data, function(x){do.call(rbind.data.frame,x)})
          duplicated_data <- do.call(rbind.data.frame,duplicated_data)
          current <- current[(!(duplicated(current[,c("id","range")]) | duplicated(current[,c("id","range")], fromLast = TRUE))),]
          current <- rbind(current,duplicated_data)
          rm(duplicated_data)
        }
        
        rm(temp)
        
        current$Region_type_numeric <- current$Region_type
        current$Region_type_numeric <- gsub("non-diploid","2",current$Region_type_numeric)
        current$Region_type_numeric <- gsub("^diploid$","1",current$Region_type_numeric)
        current$Region_type_numeric <- gsub("noisy","-1",current$Region_type_numeric)
        current$Region_type_numeric <- as.numeric(as.character(current$Region_type_numeric))
        
        cell_stats <- as.data.frame.matrix(table(data.frame(current[,c("id","Region_type")])))
        cell_stats$Cell_Ploidy <- colnames(cell_stats)[apply(cell_stats,1,which.max)]
        cell_stats$id <- row.names(cell_stats)
        current$Cell_Ploidy <- cell_stats[match(current$id, cell_stats$id),"Cell_Ploidy"]
        current$CID <- paste(current$id, current$range, sep = "-")
        temp <- current[which(tolower(current$Cell_Ploidy) == "diploid" & current$event_confidence > 100),]
        temp <- rbind(temp, current[which(tolower(current$Cell_Ploidy) != "diploid" & current$event_confidence > 50),])
        temp <- rbind(temp, current[which(current$range %in% temp$range & !current$CID %in% temp$CID),])
        temp <- reshape2::dcast(temp, formula = range~id, value.var = "Region_type_numeric")
        temp[is.na(temp)] <- 0
        row.names(temp) <- temp$range
        temp <- temp[,grep("range", colnames(temp), ignore.case = T, invert = T)]
        start_time <- Sys.time()
        temp <- as(as.matrix(temp),"CsparseMatrix")
        print("before clustering:")
        hc <- hclust(dist(t(temp)), method='complete')
        print("after clustering:")
        end_time <- Sys.time()
        print(end_time - start_time)
        hc_summary[[i]] <- hc
        print("Finnished clustering ..")
        h1names <- colnames(temp)[hc$order]
        current$id <- factor(current$id, levels = colnames(temp)[hc$order])
        cell_stats$Cell_Ploidy <- factor(cell_stats$Cell_Ploidy, levels = c("diploid","non-diploid","noisy"))
        cell_stats <- cell_stats[match(as.character(h1names),cell_stats$id),]
        cell_stats$id <- factor(cell_stats$id, levels = as.character(h1names))
        cell_stats$Mean_ploidy <- round(current_summary[match(cell_stats$id,current_summary$cell_id),"mean_ploidy"])
        current$chrom <- factor(current$chrom, levels = chr_num)
        current$Cell_type <- cell_stats[match(current$id, cell_stats$id),"Cell_Ploidy"]
        current$Cell_type <- factor(current$Cell_type, levels = c("diploid","non-diploid","noisy"))
        current$range <- factor(current$range, levels = unique(current$range))
        data_current[[i]] <- current
        names(data_current)[i] <- pheno_data[i,"SAMPLE_ID"]
        data <- rbind(data, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current))
        data_summary <- rbind(data_summary, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current_summary))
        print("data_cell_stats:")
        print(head(data_cell_stats))
        print("cell_stats:")
        print(head(cell_stats))
        data_cell_stats <- rbind(data_cell_stats, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], cell_stats))
        current <- NULL
        cell_stats <- NULL
        current_summary <- NULL
        current_mappable <- NULL
        current_group <- NULL
        print(paste("Done with ", pheno_data[i,"SAMPLE_ID"]," ..", sep = ""))
        removeModal()
      }
    }
    
    # lapply(sample_files, function(x){system(paste("rm -r ",x))})
    
    showModal(modalDialog(paste("Filtering CNV..", sep = ""), footer=NULL))
    data_group$Size <- data_group$end - data_group$start
    data_mappable$Size <- data_mappable$end - data_mappable$start
    selected_data <- data[grep("non-diploid", data$Cell_type, ignore.case = T),]
    selected_data$cell_id <- paste(selected_data$Sample, selected_data$id, sep = "_")
    selected_data_group <- data_group[which(data_group$Size >= size_threshold & data_group$event_confidence > confidence_threshold),]
    selected_data_group <- selected_data_group[which(!selected_data_group$copy_number == 2),]
    selected_data_mappable <- data_mappable[which(data_mappable$Size >= size_threshold),]
    print("dim of selected_data_mappable:")
    print(dim(selected_data_mappable))
    final_data_group <- NULL
    for(i in 1:nrow(selected_data_mappable)){
      temp <- selected_data_group[which((selected_data_group$Sample == selected_data_mappable[i,"Sample"]) &
                                          (selected_data_group$chrom == selected_data_mappable[i,"chrom"]) &
                                          (selected_data_group$start >= selected_data_mappable[i,"start"]) &
                                          (selected_data_group$end <= selected_data_mappable[i,"end"])),]
      if(nrow(temp) > 0){
        # print(i)
        final_data_group <- rbind(final_data_group, temp)
      }
      temp <- NULL
    }
    
    # filtered out CNVs in <5% of the cells
    cnv_cell_threshold <- 5/100
    total_cells <- length(unique(selected_data$cell_id))
    cnv_events <- unique(final_data_group[,grep("chrom|start|end", colnames(final_data_group), ignore.case = T)])
    final_cnv_events <- NULL
    for(i in 1:nrow(cnv_events)){
      temp <- selected_data[which((selected_data$chrom == cnv_events[i,"chrom"]) &
                                    (selected_data$start_old >= cnv_events[i,"start"]) &
                                    (selected_data$end_old <= cnv_events[i,"end"])),]
      if(nrow(temp) > 0){
        # print(i)
        temp <- temp[!duplicated(temp[,grep("cell_id|chrom|start_old|end_old", colnames(temp), ignore.case = T)]),]
        if(length(unique(temp$cell_id))/total_cells > cnv_cell_threshold){
          final_cnv_events <- rbind(final_cnv_events, cnv_events[i,])
        }
      }
      temp <- NULL
    }
    
    final_cnvs <- NULL
    print("dim of final_cnv_events:")
    print(dim(final_cnv_events))
    for(i in 1:nrow(final_cnv_events)){
      temp <- selected_data_mappable[which((selected_data_mappable$chrom == final_cnv_events[i,"chrom"]) &
                                             (selected_data_mappable$start <= final_cnv_events[i,"start"]) &
                                             (selected_data_mappable$end >= final_cnv_events[i,"end"])),]
      if(nrow(temp) > 0){
        final_cnvs <- rbind(final_cnvs, temp)
      }
      temp <- NULL
    }
    
    final_cnvs <- final_cnvs[!duplicated(final_cnvs[,grep("chrom|start|end", colnames(final_cnvs), ignore.case = T)]),grep("chrom|start|end|Size|range", colnames(final_cnvs), ignore.case = T)]
    final_cnvs$range <- paste(final_cnvs$chrom,":",final_cnvs$start,"-",final_cnvs$end, sep = "")
    
    melt_cell_binary_cnv <- NULL
    selected_cells <- unique(selected_data$cell_id)
    for(i in 1:nrow(final_cnvs)){
      temp <- selected_data[which((selected_data$chrom == final_cnvs[i,"chrom"]) &
                                    (selected_data$start_old >= final_cnvs[i,"start"]) &
                                    (selected_data$end_old <= final_cnvs[i,"end"])),]
      
      if(nrow(temp) > 0){
        temp <- unique(temp$cell_id)
        melt_cell_binary_cnv <- rbind(melt_cell_binary_cnv, data.frame(CNV_Event = final_cnvs[i,"range"], cell_id = temp , Count = 1))
      }
      temp <- NULL
    }
    
    cell_binary_cnv <- reshape2::dcast(melt_cell_binary_cnv, cell_id ~ CNV_Event, value.var = "Count")
    cell_binary_cnv[is.na(cell_binary_cnv)] <- 0
    row.names(cell_binary_cnv) <- cell_binary_cnv$cell_id
    cell_binary_cnv <- cell_binary_cnv[,grep("cell_id", colnames(cell_binary_cnv), ignore.case = T, invert = T)]
    removeModal()
    showModal(modalDialog(paste("Running dimension reduction and clustering.. ", sep = ""), footer=NULL))
    cell_groups <- find.clusters(cell_binary_cnv, max.n.clust=10, n.pca = 200, choose.n.clust = F, criterion = "goodfit")
    # 7
    # print("Fixed to 7..")
    dapc_out <- dapc(cell_binary_cnv, cell_groups$grp, n.pca = 200, n.da = 10)
    print("dapc_out:")
    print(head(dapc_out))
    dc_cnvs <- dapc_out$var.contr
    umap_out <- umap::umap(cell_binary_cnv)
    umap_coords <- data.frame(UMAP_1 = umap_out$layout[,1], UMAP_2 = umap_out$layout[,2], 
                              cell_id = row.names(umap_out$layout))
    umap_coords$Cluster <- dapc_out$assign[match(umap_coords$cell_id,row.names(dapc_out$posterior))]
    umap_coords <- umap_coords[order(umap_coords$Cluster, decreasing = F),]
    umap_coords$cell_id <- factor(umap_coords$cell_id, levels = unique(umap_coords$cell_id))
    
    hc_cnv <- hclust(dist(t(cell_binary_cnv)), method='complete')
    binary_cnv_events <- cell_binary_cnv
    binary_cnv_events$cell_id <- row.names(cell_binary_cnv)
    binary_cnv_events <- reshape2::melt(binary_cnv_events)
    colnames(binary_cnv_events) <- c("cell_id","CNV_Event","Count")
    binary_cnv_events$Count <- factor(binary_cnv_events$Count)
    binary_cnv_events$CNV_Event <- factor(binary_cnv_events$CNV_Event, levels = colnames(cell_binary_cnv)[hc_cnv$order])
    binary_cnv_events$Cluster <- dapc_out$assign[match(binary_cnv_events$cell_id,row.names(dapc_out$posterior))]
    binary_cnv_events <- binary_cnv_events[order(binary_cnv_events$Cluster, decreasing = F),]
    binary_cnv_events$cell_id <- factor(binary_cnv_events$cell_id, levels = unique(binary_cnv_events$cell_id))
    
    top_events <- data.frame(dapc_out$var.contr)
    top_events <- top_events[order(top_events$LD1, decreasing = T),]
    ld1 <- row.names(top_events)[1:50]
    top_events <- top_events[order(top_events$LD2, decreasing = T),]
    ld2 <- row.names(top_events)[1:50]
    top_events <- c(unique(ld1), unique(ld2))[1:50]
    
    top_event_chroms <- data.frame(table(gsub("(.*?):.*","\\1",top_events)))
    colnames(top_event_chroms) <- c("chrom","Count")
    top_event_chrom_threshold <- 3
    selected_chroms <- top_event_chroms[which(top_event_chroms$Count >= top_event_chrom_threshold),"chrom"]
    select_bychrom_data <- selected_data[which(selected_data$chrom %in% selected_chroms),]
    select_bychrom_data$Cluster <- dapc_out$assign[match(select_bychrom_data$cell_id,row.names(dapc_out$posterior))]
    select_bychrom_data <- select_bychrom_data[order(select_bychrom_data$Cluster, decreasing = F),]
    select_bychrom_data$CID <- paste(select_bychrom_data$cell_id, current$range, sep = "-")
    temp <- split(select_bychrom_data, select_bychrom_data$Cluster)
    for(i in 1:length(temp)){
      x <- temp[[i]]
      if(nrow(x) > 1){
        x$Ploidy <- gsub(">=10","10",x$Ploidy)
        y <- x[which(x$event_confidence > 50),]
        y <- rbind(y, x[which(x$range %in% y$range & !x$CID %in% y$CID),])
        if(nrow(y) > 100){
          x <- y
        }
        x <- dcast(x, range ~ cell_id, value.var = "Ploidy")
        row.names(x) <- x$range
        x <- x[,which(colnames(x) != "range")]
        x[is.na(x)] <- 0
        hc <- hclust(dist(t(x)), method='complete')
        x <- data.frame(Cluster = unique(temp[[i]]$Cluster), cell_id = colnames(x)[hc$order])
        temp[[i]] <- x
      }
    }
    
    temp <- do.call(rbind.data.frame, temp)
    datahcnames <- unique(temp$cell_id)
    # Workspace.RData
    select_bychrom_data <- select_bychrom_data[which(select_bychrom_data$cell_id %in% datahcnames),]
    select_bychrom_data$cell_id <- factor(select_bychrom_data$cell_id, levels = datahcnames)
    select_bychrom_data <- select_bychrom_data[,grep("cell_id|^range$|Ploidy|chrom|Cluster", colnames(select_bychrom_data), ignore.case = T)]
    select_bychrom_data <- select_bychrom_data[!duplicated(select_bychrom_data),]
    select_bychrom_data <- split(select_bychrom_data, select_bychrom_data$chrom)
    select_bychrom_data <- lapply(select_bychrom_data, function(x){
      x <- x[order(as.numeric(as.character(gsub(".*?:(.*)-(.*)","\\2",x$range)))),]
    })
    select_bychrom_data <- do.call(rbind.data.frame, select_bychrom_data)
    select_bychrom_data$range <- factor(select_bychrom_data$range, levels = unique(select_bychrom_data$range))
    # data_cell_stats$cell_id <- paste(data_cell_stats$Sample,data_cell_stats$id,sep = "_")
    # data_cell_stats$Cluster <- select_bychrom_data[match(data_cell_stats$cell_id,select_bychrom_data$cell_id),"Cluster"]
    # current <- data_cell_stats[which(!is.na(data_cell_stats$Cluster)),]
    # select_bychrom_data <- select_bychrom_data[which(select_bychrom_data$cell_id %in% data_cell_stats$cell_id),]
    # current <- current[match(unique(select_bychrom_data$cell_id), current$cell_id),]
    # current$cell_id <- factor(current$cell_id, levels = unique(current$cell_id))
    select_bychrom_data_bar <- data.frame(cell_id = select_bychrom_data$cell_id, Cluster = select_bychrom_data$Cluster)
    select_bychrom_data_bar <- select_bychrom_data_bar[match(datahcnames, select_bychrom_data_bar$cell_id),]
    select_bychrom_data_bar <- unique(select_bychrom_data_bar)
    
    data_prop <- data.frame(cell_id = unique(selected_data$cell_id), Type = "non-diploid")
    data_prop$Cluster <- dapc_out$assign[match(data_prop$cell_id,row.names(dapc_out$posterior))]
    data_prop$Sample <- gsub("(.*)_.*","\\1",data_prop$cell_id)
    data_prop <- as.data.frame.matrix(table(data_prop[,c("Sample","Cluster")]))
    current <- data[grep("noisy", data$Cell_type, ignore.case = T, invert = T),c("Sample","id")]
    current <- current[!duplicated(current),]
    total_cells <- data.frame(table(current[,c("Sample")]))
    total_cells <- total_cells[match(row.names(data_prop), total_cells$Var1),]
    for(i in 1:nrow(data_prop)){
      data_prop[i,] <- data_prop[i,]/total_cells[i,"Freq"]
    }
    data_prop$Sample <- row.names(data_prop)
    data_prop <- melt(data_prop)
    colnames(data_prop) <- c("Sample","Cluster","Proportion")
    
    cn_clusters <- data.frame(cell_id = unique(selected_data$cell_id))
    cn_clusters$Cluster <- dapc_out$assign[match(cn_clusters$cell_id,row.names(dapc_out$posterior))]
    cn_clusters$copy_number <- round(data_summary$mean_ploidy)[match(cn_clusters$cell_id, paste(data_summary$Sample,data_summary$cell_id,sep = "_"))]
    cn_clusters <- cn_clusters[which(!is.na(cn_clusters$Cluster) & !is.na(cn_clusters$copy_number)),]
    dcast_cn_clusters <- reshape2::dcast(cn_clusters, cell_id ~ Cluster, value.var = "copy_number")
    dcast_cn_clusters[is.na(dcast_cn_clusters)] <- 0
    row.names(dcast_cn_clusters) <- dcast_cn_clusters$cell_id
    dcast_cn_clusters <- dcast_cn_clusters[,grep("cell_id", colnames(dcast_cn_clusters), ignore.case = T, invert = T)]
    tree_est <- fastme.bal(dist(t(dcast_cn_clusters)))
    
    plotx <- split(cn_clusters, cn_clusters$Cluster)
    plotx <- lapply(plotx, function(x){
      x <- data.frame(Cluster = unique(x$Cluster),
                      Median_CN = median(x$copy_number),
                      Size = nrow(x))
    })
    plotx <- do.call(rbind.data.frame,plotx)
    plotx$Cluster <- as.character(plotx$Cluster)
    plotx <- plotx[match(tree_est$tip.label, plotx$Cluster),]
    clusters <- split(plotx[,"Cluster"], plotx$Cluster)
    tree_est <- groupOTU(tree_est, clusters)
    tree_est$plotx <- plotx
    
    results <- NULL
    results$data_current <- data_current
    results$data_cell_stats <- data_cell_stats
    print("saving dapc_out:")
    print(head(dapc_out))
    results$dapc_out <- dapc_out
    results$umap_coords <- umap_coords
    results$binary_cnv_events <- binary_cnv_events
    results$select_bychrom_data <- select_bychrom_data
    results$select_bychrom_data_bar <- select_bychrom_data_bar
    results$data_prop <- data_prop
    results$tree_est <- tree_est
    results$sample_colors <- gen_colors(color_conditions$tenx, length(unique(data_cell_stats$Sample)))
    names(results$sample_colors) <- unique(data_cell_stats$Sample)
    results$cluster_colors <- gen_colors(color_conditions$bright, length(unique(umap_coords$Cluster)))
    names(results$cluster_colors) <- sort(as.numeric(as.character(unique(umap_coords$Cluster))))
    results$project_name <- project_name
    system(paste("rm -r ",cdir, sep = ""))
    ##################################################################################################
    
    # saveRDS(results,"results.RDS")
    annot_names <- names(results$data_current)
    annot_names <- sort(annot_names)
    removeModal()
    
    updateSelectInput(session, inputId = 'p1id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    print("Completed!")
    return(results)
  })
  
  output$p1plot <- renderPlot({ #renderPlotly
    if(input$p1id != ""){
      showModal(modalDialog("Plotting figure 1..", footer=NULL))
      current <- data()[['data_current']][[which(names(data()[['data_current']]) == input$p1id)]]
      cell_stats <- data()[['data_cell_stats']][which(data()[['data_cell_stats']]$Sample == input$p1id),]
      cell_stats$id <- factor(cell_stats$id, levels = levels(current$id))
      print("Printing ploidy by chrom position plot ..")
      p1 <- ggplot(current, aes(range, id, color = Ploidy, fill = Ploidy))+ # [sample(1:nrow(data),size = 100, replace = F),]
        theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),axis.text.y=element_blank(),
              plot.title = element_text(size=25, face = "bold")) +
        # legend.text=element_text(size=18)) + # hjust = 0.5
        # geom_point(alpha = 0.6)+
        geom_tile(size = 0.8)+
        xlab("Chromosome Positions") + ylab("") +
        scale_color_manual(values = color_heatmap)+
        scale_fill_manual(values = color_heatmap)+
        facet_grid(.~chrom, scales = "free", switch = "x", space = "free_x")+
        ggtitle(paste(input$p1id, ", ", signif(length(which(current$Ploidy == "2"))/nrow(current)*100,3), "% diploid"))
      
      p2 <- ggplot(cell_stats)+
        geom_bar(mapping = aes(x = 1, y = id, fill = Cell_Ploidy), stat = "identity", width = 1)+
        theme_void()+
        theme(panel.spacing.x = unit(1, "mm")) +
        scale_fill_manual(values = cell_ploidy_colors)
      
      legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
      p1 <- p1 + theme(legend.position = "none")
      p2 <- p2 + theme(legend.position = "none")
      p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))
      p <- plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1))
      removeModal()
      return(p)
    }
  }, height = 500, width = 900)
  
  observeEvent(input$p1plot, {
    screenshot(id="p1plot", filename = paste("1SCA_PLOIDY_INFO_CHROM_POSITION_",input$p1id, sep = ""), scale = 2)
  })
  
  output$p2plot <- renderPlot({
    showModal(modalDialog("Plotting figure 2..", footer=NULL))
    plotx <- data()[['data_cell_stats']]
    p <- ggplot(plotx, aes(x = Mean_ploidy, y = Sample)) +
      geom_density_ridges(aes(fill = Sample)) + ggtitle(data()$project_name) +
      xlab("Mean single-cell ploidy") + ylab("Samples") +
      scale_fill_manual(values = data()[['sample_colors']])
    p <- adjust_theme(p, xsize = 20, title_size = 25)
    removeModal()
    return(p)
  }, height = 800, width = 1200)
  
  observeEvent(input$p2plot, {
    screenshot(id="p2plot", filename = paste("2SCA_PLOIDY_INFO_MEAN_CELL_PLOIDY_CNV_",input$p2id, sep = ""), scale = 2)
  })
  
  output$p3plot <- renderPlot({
    showModal(modalDialog("Plotting figure 3..", footer=NULL))
    plotx <- data()[['dapc_out']]
    p <- scatter(plotx, bg="white", col = data()[['cluster_colors']], legend=T,
                 scree.da=FALSE, inset.solid = 0.6,cex.lab = 1, label.inds = c("DAPC_1","DAPC_2"))
    removeModal()
    return(p)
  }, height = 800, width = 1200)
  
  observeEvent(input$p3plot, {
    screenshot(id="p3plot", filename = paste("3SCA_DAPC_COMPONENT12_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p4plot <- renderPlot({
    showModal(modalDialog("Plotting figure 4..", footer=NULL))
    plotx <- data()[['umap_coords']]
    p <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "Cluster", plot_title = data()[['project_name']],
                      col = data()[['cluster_colors']], annot = T, legend_position = "right", numeric = T, 
                      point_size = 2, label_size = 8,legendsize = 15)
    removeModal()
    return(p)
  }, height = 800, width = 1200)
  
  observeEvent(input$p4plot, {
    screenshot(id="p4plot", filename = paste("4SCA_UMAP_COMPONENT12_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p5plot <- renderPlot({
    showModal(modalDialog("Plotting figure 5..", footer=NULL))
    plotx <- data()[['binary_cnv_events']]
    p1 <- ggplot(plotx, aes(CNV_Event, cell_id, fill = Count))+
      theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
            axis.ticks.y=element_blank(),axis.text.y=element_blank(),
            plot.title = element_text(size=25, face = "bold"),
            axis.title = element_text(size=20, face = "bold"),
            legend.title = element_text(size =20, face = "bold"),
            legend.text = element_text(size = 15)) +
      geom_tile()+
      xlab("CNV Events") + ylab("") +
      scale_fill_manual(values = c("black","#eec643"))
    
    p2 <- ggplot(plotx)+
      geom_bar(mapping = aes(x = 1, y = cell_id, fill = Cluster), stat = "identity", width = 1)+
      theme_void()+
      theme(panel.spacing.x = unit(1, "mm"),
            legend.title = element_text(size =20, face = "bold"),
            legend.text = element_text(size = 15)) +
      scale_fill_manual(values = data()[['cluster_colors']])
    
    legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
    p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))
    p <- plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1))
    removeModal()
    return(p)
  }, height = 800, width = 1200)
  
  observeEvent(input$p5plot, {
    screenshot(id="p5plot", filename = paste("5SCA_BINARY_CNV_EVENTS_BY_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p6plot <- renderPlot({
    showModal(modalDialog("Plotting figure 6..", footer=NULL))
    plotx <- data()[['select_bychrom_data']]
    plotx_bar <- data()[['select_bychrom_data_bar']]
    
    p1 <- ggplot(plotx, aes(range, cell_id, color = Ploidy, fill = Ploidy))+
      theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
            axis.ticks.y=element_blank(),axis.text.y=element_blank(),
            plot.title = element_text(size=25, face = "bold"),
            axis.title = element_text(size=20, face = "bold"),
            legend.title = element_text(size =20, face = "bold"),
            legend.text = element_text(size = 15)) +
      geom_tile(size = 0.8)+
      scale_fill_manual(values = color_heatmap)+
      scale_color_manual(values = color_heatmap)+
      xlab("Chromosome Positions") + ylab("") +
      facet_grid(.~chrom, scales = "free", switch = "x", space = "free_x")
    
    p2 <- ggplot(plotx_bar)+
      geom_bar(mapping = aes(x = 1, y = cell_id, fill = Cluster), stat = "identity", width = 1)+
      theme_void()+
      theme(panel.spacing.x = unit(1, "mm"),
            legend.title = element_text(size =20, face = "bold"),
            legend.text = element_text(size = 15)) +
      scale_fill_manual(values = data()[['cluster_colors']])
    
    legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
    p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))
    p <-  plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1))
    removeModal()
    return(p)
  }, height = 800, width = 1400)
  
  observeEvent(input$p6plot, {
    screenshot(id="p6plot", filename = paste("6SCA_PLOIDY_INFO_TOP_50_CNVs_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p7plot <- renderPlot({
    showModal(modalDialog("Plotting figure 7..", footer=NULL))
    plotx <- data()[['data_prop']]
    p <- ggplot(plotx, aes(x=Sample, y=Proportion, fill=Cluster, group = Cluster)) + 
      geom_area(alpha=0.9, size=0.5, colour = "black")+
      theme_classic()+ggtitle(data()[['project_name']])+
      ylab("Fraction of Cells")+
      scale_fill_manual(values = data()[['cluster_colors']])+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p <- adjust_theme(p, legend = "right", xsize = 20, xangle = 45, hejust = 1, vejust = 1)
    removeModal()
    return(p)
    
  }, height = 800, width = 1400)
  
  observeEvent(input$p7plot, {
    screenshot(id="p7plot", filename = paste("7SCA_CLUSTER_PROPORTION_BY_SAMPLES_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p8plot <- renderPlot({
    showModal(modalDialog("Plotting figure 8..", footer=NULL))
    tree_est <- data()[['tree_est']]
    p <- ggtree(tree_est, layout="daylight",
                branch.length = c(tree_est$plotx$Size/300,rep(min(tree_est$plotx$Size/300),(nrow(tree_est$edge)+1) - nrow(tree_est$plotx))),
                aes(color = group,
                    size = c(tree_est$plotx$Size/50,rep(min(tree_est$plotx$Size/100),(nrow(tree_est$edge)+1) - nrow(tree_est$plotx))))) + 
      geom_tiplab(hjust = -2, offset=.1) +
      scale_color_manual(values = data()[['cluster_colors']]) +
      theme(legend.position="none") +
      scale_size_continuous()
    removeModal()
    return(p)
    
  }, height = 1200, width = 1400)
  
  observeEvent(input$p8plot, {
    screenshot(id="p8plot", filename = paste("8SCA_PHYLOGENETIC_TREE_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
  })
})
