#!/usr/bin/env Rscript
##############################################################################
#             SCA CyTOF Analysis Package
#             SingleCellAnalyst.org Production
#             Purpose: Multi-Omics Single Cell Analysis Pipeline Functions
#             Version: V1.1.0
#             Author: Lu Pan
#             Last-update Date: 2023-05-15
##############################################################################

Packages <- c("Seurat","ggplot2","patchwork","dplyr","ggridges",
              "ggthemes","gridExtra","reshape2", "ggrepel","scales","xlsx")

suppressMessages(lapply(Packages, library, character.only = TRUE))

gen_colors <- function(col = NULL,n, lum = "light"){
  if(is.null(col)){
    col <- circlize::rand_color(n = n, luminosity = lum)
  }else{
    col <- colorRampPalette(col)(n)
  }
  return(col)
}

show_colors <- function(cols){
  return(show_col(cols))
}
color_ini <- function(){
  color_conditions <- NULL
  # color_conditions$ggplot <- gg_color_hue(20)
  color_conditions$monet <- c("#76A2BB","#BCB180","#E8AC75","#EE9873","#D69896","#F4AFA3","#1166AE","#1C9FD3","#C45A35","#9C9B0B","#E34615","#F6E576","#1A3D40","#7F2125","#1E6652","#1BA0C7","#075B25","#0F8441","#F6C04A","#E1B5B2")
  color_conditions$vanngogh <- c("#EEE6B2","#F9E956","#A66100","#54ABAB","#4062B4","#6D87B9","#E6AE95","#5D1C01","#964000","#010101","#F15E01","#EB9D00","#5B1700","#6E1700","#8A7B01","#DED6AD")
  color_conditions$manycolors <- c(unique(unlist(lapply(ggthemes::ggthemes_data$tableau$`color-palettes`$regular,function(x){x <- x$value}))))
  color_conditions$ggplot <- c("#F8766D", "#EA8331", "#D89000", "#C09B00", "#A3A500", "#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3", "#00BFC4", "#00BAE0", "#00B0F6", "#35A2FF", "#9590FF", "#C77CFF", "#E76BF3", "#FA62DB", "#FF62BC", "#FF6A98")
  color_conditions$colorful <- c("#EC6077", "#9AF270", "#CBD157", "#F0914C", "#6745A4", "#70E695", "#56B6CD", "#5076D7", "#B94AA7")
  color_conditions$general <- c("#5F75A5","#91D376","#D75B58","#F5BFD3","#A8C9EA","#B09C16","#F69F99","#AC79A3","#E89314","#EAD256",
                                "#78706E","#D1A5CA","#F7C277","#569794","#B9B0AC","#99785E","#5FA346","#8DBCB6","#CC7296","#D3B6A5")
  color_conditions$tableau20 <- c("#9edae5","#17becf","#dbdb8d","#bcbd22","#c7c7c7","#7f7f7f","#f7b6d2","#e377c2","#c49c94","#8c564b",
                                  "#c5b0d5","#9467bd","#ff9896","#d62728","#98df8a","#2ca02c","#ffbb78","#ff7f0e","#aec7e8","#1f77b4")
  color_conditions$tenx <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold", "#a65628", "#999999", "#9e2a2b", "grey", "#58bc82", "purple")
  color_conditions$mark <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                             "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                             "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                             "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                             "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                             "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d")
  color_conditions$bright <- c("#6C85B6","#F5B317","#E57369","#94C1CB","#61AE87","#BFC029","#B67694")
  color_conditions$cold <- c("#5F75A5","#91d376","#A8C9EA","#78706E","#569794","#8DBCB6","#5FA346","#B9B0AC")
  color_conditions$warm <- c("#D75B58","#E89314","#F5BFD3","#F69F99","#F7C277","#CC7296","#D3B6A5","#D1A5CA","#99785E","#B09C16","#EAD256","#AC79A3")
  color_conditions$alternate <- c("#5F75A5","#E89314","#91d376","#D75B58","#A8C9EA","#F5BFD3","#78706E","#F69F99","#569794","#F7C277","#8DBCB6","#CC7296","#5FA346","#D3B6A5","#B9B0AC","#D1A5CA")
  color_conditions$gradient <- c("#584C8E","#5079A9","#5EA1AC","#77C1A1","#A2CF9F",
                                 "#CBDE9C","#E9ECA5","#FAE79E","#F7CC82","#EEA56C",
                                 "#E77A58","#D4534E","#B2334A","#8A2140")
  color_conditions$redbluecont <- c("#26456e","#1c5998","#1c73b1","#3a87b7","#67add4","#cacaca","#fc8375","#df513f","#d11719","#bd1316","#9c0824")
  color_conditions$orangebluecont <- c("#26456e","#1c5998","#1c73b1","#3a87b7","#67add4","#7bc8e2","#cacaca","#fdab67","#fd8938","#f06511","#d74401","#a33202","#7b3014")
  color_conditions$trafficlightcont <- c("#9fcd99","#ffdd71","#f26c64")
  color_conditions$orangegradient <- c("white", "cornflowerblue", "yellow", "red")
  color_conditions$flow <- c("blue","cyan","green","yellow","orange","red","red4")
  color_conditions$Rainbow <- rev(rainbow(10))
  color_conditions$OrangeBlue <- c("#26456e","#1c5998","#1c73b1","#3a87b7","#67add4","#7bc8e2","#cacaca","#fdab67", "#fd8938","#f06511","#d74401","#a33202","#7b3014")
  color_conditions$GreenRed <- c("green", "white", "red")
  color_conditions$BlueRed <- c("blue", "#EEEEEE", "red")
  color_conditions$Viridis <- hcl.colors(10)
  color_conditions$Heat <- heat.colors(10)
  # color_conditions$Roma <- hcl.colors(10, palette = "Roma")
  color_conditions$Plasma <- hcl.colors(10, palette = "Plasma")
  color_conditions$Zissou <- hcl.colors(10, palette = "Zissou 1")
  color_conditions$RedYellowBlue <- hcl.colors(10, palette = "RdYlBu")
  color_conditions$Spectral <- hcl.colors(10, palette = "Spectral")
  color_conditions$Terrain <- terrain.colors(10)
  color_conditions$CM <- cm.colors(10)
  color_conditions$BlueYellowRed <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
  color_conditions$solarized <- ggthemes::ggthemes_data$solarized
  color_conditions$ptol <- ggthemes::ggthemes_data$ptol
  color_conditions$tableau <- ggthemes::ggthemes_data$tableau$`color-palettes`$regular
  color_conditions$tableaudivergent <- ggthemes::ggthemes_data$tableau$`color-palettes`$`ordered-diverging`
  color_conditions$tableausequential <- ggthemes::ggthemes_data$tableau$`color-palettes`$`ordered-sequential`
  color_conditions$ggthemeothers <- ggthemes::ggthemes_data$others
  
  
  return(color_conditions)
}

shiny_cols <- function(cols,n){
  color_conditions <- color_ini()
  color_conditions <- color_conditions[1:9]
  names(color_conditions) <- c("Default","Colorful","Tableau1","Tableau2","Tenx","Mark","Bright","Cold","Warm")
  col_out <- gen_colors(col = color_conditions[[cols]], n)
  return(col_out)
}

readfile <- function(file_dir){
  library("xlsx")
  
  if(length(grep("\\.csv$", file_dir, ignore.case = T)) > 0){
    if(ncol(read.csv(file_dir, header = T, sep = ";")) > 1){
      current <- suppressWarnings(read.csv(file_dir, header = T, sep = ";"))
    }else{
      current <- suppressWarnings(read.csv(file_dir, header = T))
    }
  }else if(length(grep("\\.txt$", file_dir, ignore.case = T)) > 0){
    if(ncol(suppressWarnings(read.table(file_dir, header = T))) == 1){
      current <- suppressWarnings(read.table(file_dir, header = T, sep = ";"))
    } else if(ncol(suppressWarnings(read.table(file_dir, header = T, sep = ","))) > 1){
      current <- suppressWarnings(read.table(file_dir, header = T, sep = ","))
    }else{
      current <- suppressWarnings(read.table(file_dir, header = T))
    }
  }else if(length(grep("\\.xls$|\\.xlsx$", file_dir, ignore.case = T)) > 0){
    current <- suppressWarnings(read.xlsx(file_dir, 1, header=TRUE))
  }else{
    current <- suppressWarnings(read.table(file_dir, header=TRUE, sep = "\t"))
  }
  
  if(toupper(colnames(current)[1]) == "X"){
    row.names(current) <- current[,1]
    current <- current[,which(toupper(colnames(current)) != "X")]
  }
  return(current)
}

pheno_ini <- function(phenodata_dir, pipeline, extra_para = NULL, isDir = TRUE){
  
  if(isDir == TRUE){
    pheno_data <- readfile(phenodata_dir)
  }else{
    pheno_data <- phenodata_dir
  }
  if(toupper(pipeline) == "SPATIAL"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","GROUP")
  }else if (toupper(pipeline) == "SCRNASEQ"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP") # "CELL_CYCLE", "DATA_TYPE","CELL_TYPE","PAIR_ID","TISSUE_TYPE","CANCER_TYPE"
  } else if(toupper(pipeline) == "SMARTSEQ2"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
  }else if(toupper(pipeline) == "CNV"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID")
  }else if(toupper(pipeline) == toupper("scIMMUNE")){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP","DATA_TYPE","CELL_TYPE","PAIR_ID")
  }else if(toupper(pipeline) == "BULK"){
    if(length(grep("Treatment", colnames(pheno_data), ignore.case = TRUE)) > 0){
      colnames(pheno_data) <- c("SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP","TREATMENT")
    }else{
      colnames(pheno_data) <- c("SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
    }
  }else if(toupper(pipeline) == "BULK_TCR"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","CELL_TYPE","GROUP","BATCH")
  }else if (toupper(pipeline) == "FLOW"){
    # if(!is.null(extra_para)){
    #   colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP",extra_para)
    # }else{
      colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP")
    # }
  }else if (toupper(pipeline) == "ATAC"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP","REF_GENOME","RNASEQ_FILE","DATA_TYPE")
  }else if (toupper(pipeline) %in% c("EV","CYTOF")){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
    # if(!is.null(extra_para)){
    #   colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
    # }else{
    #   colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
    # }
  }
  if((length(grep("run|BATCH", pheno_data$BATCH, ignore.case = T)) == 0) &
     (length(grep("BATCH", colnames(pheno_data), ignore.case = T)) > 0) &
     (all(grepl("^[0-9]+$", pheno_data$BATCH)))){
    pheno_data$BATCH <- paste("BATCH_", pheno_data$BATCH, sep = "")
  }
  
  if((length(grep("SAMPLE|^S_|^S.*", pheno_data$SAMPLE_ID, ignore.case = T)) == 0) &
     (length(grep("SAMPLE_ID", colnames(pheno_data), ignore.case = T)) > 0) &
     (all(grepl("^[0-9]+$", pheno_data$SAMPLE_ID)))){
    pheno_data$SAMPLE_ID <- paste("SAMPLE_", pheno_data$SAMPLE_ID, sep = "")
  }
  
  if((length(grep("INDIVIDUAL_ID", colnames(pheno_data), ignore.case = T)) > 0) &
     (all(grepl("^[0-9]+$", pheno_data$INDIVIDUAL_ID)))){
    pheno_data$INDIVIDUAL_ID <- paste("P", pheno_data$INDIVIDUAL_ID, sep = "")
  }
  
  for(i in grep("file|ref_genome|RNASEQ_FILE", colnames(pheno_data), ignore.case = T, invert = T)){
    pheno_data[,i] <- gsub("^\\s+|\\s+$|\\.fcs|\\.csv|\\.h5|\\.txt|\\.mtx|\\.tsv|\\.gz","",pheno_data[,i], ignore.case = T)
    pheno_data[,i] <- gsub("\\s+","_",pheno_data[,i], ignore.case = T)
    pheno_data[,i] <- gsub("-","_",pheno_data[,i], ignore.case = T)
    pheno_data[,i] <- toupper(pheno_data[,i])
  }
  
  if(toupper(pipeline) == "BULK"){
    for(i in grep("SAMPLE_ID", colnames(pheno_data), ignore.case = T, invert = T)){
      pheno_data[,i] <- factor(pheno_data[,i])
    }
  }
  
  # if((length(grep("TISSUE_TYPE|CANCER_TYPE", colnames(pheno_data), ignore.case = T)) > 0)){
  #   for(i in grep("TISSUE_TYPE|CANCER_TYPE", colnames(pheno_data), ignore.case = T)){
  #     pheno_data[,i] <- capFirstletter(pheno_data[,i])
  #     if(length(which(is.na(pheno_data[,i]))) > 0){
  #       pheno_data[which(is.na(pheno_data[,i])),i] <- "NULL"
  #     }
  #     if(length(grep("NULL",pheno_data[,i], ignore.case = T))>0){
  #       pheno_data[grepl("NULL",pheno_data[,i], ignore.case = T) | 
  #                    is.null(pheno_data[,i]) | (pheno_data[,i] == "") |
  #                    toupper(pheno_data[,i]) == "NA",i] <- "NULL"
  #     }
  #     pheno_data[,i] <- gsub("\\.|-|_"," ",pheno_data[,i])
  #     if(length(grep("NO.*ANNOTATION",pheno_data[,i], ignore.case = T))>0){
  #       pheno_data[grep("NO.*ANNOTATION",pheno_data[,i], ignore.case = T),i] <- "NO_ANNOTATION"
  #     }
  #   }
  # }
  
  return(pheno_data) 
}



capFirstletter <- function(line){
  line <- paste(toupper(substring(line,1,1)),tolower(substring(line,2,nchar(line))), sep = "")
  return(line)
}

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

prepare_cols <- function(col, n){
  
  color_conditions <- color_ini()
  if(length(col) > 1 & !is.null(col)){
    col <- colorRampPalette(col)(n)
  } else if(is.null(col)){
    col <- colorRampPalette(color_conditions$colorful)(n)
  }else if(toupper(col) == toupper("default")){
    col <- gg_color_hue(n)
  }
  
}

add_names <- function(current, sample_name,current_ident){
  sample_name <- gsub("_FILTERED_PEAK_BC|_FILTERED_FEATURE_BC","", sample_name, ignore.case = T)
  # sample_name <- ifelse(nchar(sample_name) > 8, gsub(".*\\.(.*)","\\1",sample_name), sample_name)
  current <- SetIdent(current, value = current_ident)
  levels(current@active.ident) <- sample_name
  current@project.name <- sample_name
  current$orig.ident <- sample_name
  return(current)
}

adjust_theme <- function(p, xangle = 0,legend = "right", title_size = 20, xsize=20, hejust = 0, vejust = 0, strip_size = 20){
  p <- p+ theme_classic()+
    # ggtitle(title) +
    # scale_fill_manual(values = col) +
    theme(axis.text.x = element_text(size = xsize, angle = xangle, hjust=hejust,vjust = vejust),
          axis.text.y = element_text(size = 20), legend.position = legend,
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = strip_size),
          strip.background = element_blank(),
          plot.title = element_text(size =title_size, face = "bold", hjust = 0.5))
  return(p)
}

create_centroids <- function(plotx, x, y, group){
  
  centroids <- split(plotx, plotx[,group])
  # centroid_groups <- names(centroids)
  centroids <- lapply(centroids, function(cent){
    cent <- data.frame(Cluster = ifelse(nrow(cent) > 0,as.character(unique(cent[,group])), NA),
                       Centroid_X = ifelse(nrow(cent) > 0, median(cent[,x]), NA),
                       Centroid_Y = ifelse(nrow(cent) > 0, median(cent[,y]), NA))})
  centroids <- do.call(rbind.data.frame, centroids)
  # centroids[,group] <- centroid_groups
  centroids <- centroids[!is.na(centroids[,"Cluster"]),]
  return(centroids)
  
}

plot_bygroup <- function(plotx, x, y, group, plot_title, col = NULL, annot = TRUE, legend_position = "right", numeric = FALSE,point_size = 4, label_size = 5, legendsize = 20){
  
  color_conditions <- color_ini()
  n <- length(unique(plotx[,group]))
  if(is.null(names(col)) | n > length(col)){
    # if(n > 1){
    col <- prepare_cols(col, n)
    # }
  }else{
    if(!is.null(names(col))){
      col <- col[which(names(col) %in% unique(plotx[,group]))]
    }else{
      col <- col[1:length(unique(plotx[,group]))]
    }
  }
  
  if(annot == TRUE){
    if(numeric == FALSE){
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.character(plotx[,group]))))
    }else{
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.numeric(as.character(plotx[,group])))))
    }
    centroids <- create_centroids(plotx, x, y, group)
    centroids$Col <- col
    
    plotx$Label <- ""
    for(i in 1:nrow(centroids)){
      plotx[nrow(plotx)+1,x] <- centroids[i,"Centroid_X"]
      plotx[nrow(plotx),y] <- centroids[i,"Centroid_Y"]
      plotx[nrow(plotx),group] <- centroids[i,"Cluster"]
      plotx[nrow(plotx),"Label"] <- centroids[i,"Cluster"]
    }
    
    p <- ggplot(plotx, aes(x = plotx[,x], y = plotx[,y], color = plotx[,group], label = Label)) + 
      geom_point(alpha = ifelse(plotx$Label != "", 0, 1), size = point_size) +
      scale_color_manual(values = col) +
      theme_classic() +
      xlab(x) + ylab(y) + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 25),
            legend.position = legend_position,
            legend.title = element_blank(),
            legend.text = element_text(size = legendsize),
            legend.key.size = unit(0.8, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(color=guide_legend(override.aes = list(size = 3), title=toupper(group), ncol=ifelse(length(unique(plotx[,group])) > 20, 2, 1)))
    
    if(numeric == FALSE){
      p <- p + geom_text_repel(max.overlaps = Inf,force=1,
                               point.padding = 0, # additional padding around each point
                               min.segment.length = 0, # draw all line segments
                               max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
                               box.padding = 0.3, size = label_size, colour = "black")
      # p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = label_size, hjust = 0, fontface =1)
    }else{
      p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = label_size, hjust = 0, fontface =1)
    }
  }else{
    p <- ggplot(plotx, aes(x = plotx[,x], y = plotx[,y], color = plotx[,group])) + 
      geom_point(alpha = 1, size = point_size) +
      scale_color_manual(values = col) +
      theme_classic() +
      xlab(x) + ylab(y) + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 25),
            legend.position = legend_position,
            legend.title = element_blank(),
            legend.text = element_text(size = legendsize),
            legend.key.size = unit(1, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(color=guide_legend(override.aes = list(size = 3), title=toupper(group), ncol=ifelse(length(unique(plotx[,group])) > 20, 2, 1)))
  }
  
  return(p)
  
}

gen10x_plotx <- function(data, groups = NULL, selected = "ALL", include_meta = FALSE){
  if(toupper(selected[1]) == "ALL"){
    plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                        tSNE_1 = data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                        tSNE_2 = data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                        PC_1 = data@reductions$pca@cell.embeddings[,"PC_1"],
                        PC_2 = data@reductions$pca@cell.embeddings[,"PC_2"])
  }else{
    for(i in 1:length(selected)){
      if(i == 1){
        if(toupper(selected[i]) == "UMAP_SELECTED"){
          plotx <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("_SELECTED","",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }else{
          plotx <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("PCA|PCA_SELECTED","PC",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }
        colnames(plotx) <- paste(toupper(gsub("PCA","PC",selected[i], ignore.case = T)),c("_1","_2"), sep = "")
      }else{
        if(toupper(selected[i]) == "UMAP_SELECTED"){
          temp <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("_SELECTED","",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }else{
          temp <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("PCA|PCA_SELECTED","PC",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }
        colnames(temp) <- paste(toupper(gsub("PCA","PC",selected[i], ignore.case = T)),c("_1","_2"), sep = "")
        plotx <- cbind(plotx, temp)
      }
    }
    colnames(plotx) <- toupper(gsub("PCA_([0-9]+)$","PC_\\1",colnames(plotx), ignore.case = T))
    colnames(plotx) <- toupper(gsub("tSNE","tSNE",colnames(plotx), ignore.case = T))
  }
  
  if(!is.null(groups)){
    for(i in 1:length(groups)){
      plotx[,groups[i]] <- data@meta.data[,groups[i]]
    }
  }
  
  if(include_meta == TRUE){
    plotx <- cbind(data@meta.data, plotx)
  }
  
  return(plotx)
  
}

showCols <- function(cl=colors(), bg = "grey", cex = 0.75, rot = 30){
  # showCols(bg="gray20",cl=brewer.pal(11, "RdYlBu"), rot=30, cex=0.9)
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}

scale_plotx <- function(plotx, direction = "BOTH"){
  plotx <- plotx[,grep("Sample|EV",colnames(plotx), ignore.case = T, invert = T)]
  if(toupper(direction) == "ROW"){
    for(i in 1:nrow(plotx)){
      plotx[i,] <- plotx[i,]/sum(plotx[i,])
    }
  } else if (toupper(direction) == "COLUMN"){
    for(i in 1:ncol(plotx)){
      plotx[,i] <- plotx[,i]/sum(plotx[,i])
    }
  } else if (toupper(direction) == "BOTH"){
    for(i in 1:nrow(plotx)){
      plotx[i,] <- plotx[i,]/sum(plotx[i,])
    }
    
    for(i in 1:ncol(plotx)){
      plotx[,i] <- plotx[,i]/sum(plotx[,i])
    }
  }
  plotx <- log2(plotx+1)
  plotx <- scale(plotx)
  plotx <- t(scale(t(plotx)))
  plotx <- plotx[,grep("NA",colSums(plotx),ignore.case = T,invert = T)]
  plotx <- plotx[,which(!colSums(plotx) == 0)]
  return(plotx)
}

sca_ridge <- function(plotx, x, y, fill, cols){
  p <- ggplot(plotx, aes_string(x = x, y = y)) +
    geom_density_ridges(aes_string(fill = fill)) +
    xlab(x) + ylab(y) +
    scale_fill_manual(values = shiny_cols(cols, length(unique(plotx[,fill]))))+
    theme_classic()
  p <- adjust_theme(p)
  return(p)
}

flow_plot <- function(data, chnl1,chnl2, bin = 300){
  p <- autoplot(data, chnl1,chnl2, bin = 350, strip.text = "gate") + 
    geom_gate(g, colour = "red", size = 0.7) + 
    geom_stats(adjust = 0.5) + theme_classic() + 
    scale_fill_gradientn(colours = c("blue","cyan","green","yellow","orange","red","red4"))+
    theme(strip.text = element_text(size = 20, face = "bold"), 
          strip.background = element_blank(),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  return(p)
}

complex_heatmap <- function(x, legend_title = NULL, col_title = NULL, row_title = NULL, 
                            col = rev(rainbow(10)),nrow_clus = 1, row_annot = NULL, legendtitle = "EXPRESSION"){
  library("ComplexHeatmap")
  p <- Heatmap(x, name = legend_title, col = col, row_km = nrow_clus, right_annotation = row_annot,
               row_names_gp = gpar(fontsize = 8), row_names_side = "left",
               column_title = col_title,column_names_rot = 45,
               heatmap_legend_param = list(title = legendtitle))
  return(p)
}

own_facet_scatter <- function(plotx, feature1, feature2, isfacet = T, 
                              title, col=color_conditions$bright, color_by,
                              group_by = NULL, xlabel = feature1, ylabel = feature2,
                              strip_size = 15, legend_pos = "right", ncol = 2){
  # set.seed(2022)
  p <- ggplot(plotx, aes(x = plotx[,feature1], y = plotx[,feature2], color = plotx[,color_by])) + 
    geom_point(size = 0.9) +
    scale_color_manual(values = gen_colors(col, length(unique(plotx[,color_by])))) +
    # scale_color_manual(values = sample(gen_colors(col, length(unique(plotx[,color_by]))), size = length(unique(plotx[,color_by])), replace = F)) +
    theme_classic() + 
    xlab(xlabel)+ylab(ylabel)+
    ggtitle(title) +
    guides(colour = guide_legend(title = color_by, override.aes = list(size=5)))+
    theme(legend.position = legend_pos, 
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.8, "cm"),
          plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
          strip.text = element_text(size = strip_size, face = "bold"),
          strip.background = element_blank())
  if(isfacet == T & !is.null(group_by)){
    p <- p + facet_wrap(~plotx[,group_by], ncol = ncol)
  }
  return(p)
}
