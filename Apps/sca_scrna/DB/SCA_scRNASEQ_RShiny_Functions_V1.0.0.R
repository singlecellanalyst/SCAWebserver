#!/usr/bin/env Rscript
##############################################################################
#             SCA scCNV Analysis Package
#             SingleCellAnalyst.org Production
#             Purpose: Multi-Omics Single Cell Analysis Pipeline Functions
#             Version: V1.0.0
#             Author: Lu Pan
#             Date: 2020-02-01
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

color_conditions <- color_ini()

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
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP") # "CELL_CYCLE", "DATA_TYPE","CELL_TYPE","PAIR_ID","TISSUE_TYPE","CANCER_TYPE"
  } else if(toupper(pipeline) == "SMARTSEQ2"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
  }else if(toupper(pipeline) == "CNV"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID")
  }else if(toupper(pipeline) == toupper("scIMMUNE")){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP","DATA_TYPE","CELL_TYPE","PAIR_ID")
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

genharmony_plotx <- function(data, groups = NULL){
  plotx <- data.frame(PC_1 = data@reductions$pca@cell.embeddings[,"PC_1"],
                      PC_2 = data@reductions$pca@cell.embeddings[,"PC_2"],
                      HARMONY_1 = data@reductions$harmony@cell.embeddings[,"harmony_1"],
                      HARMONY_2 = data@reductions$harmony@cell.embeddings[,"harmony_2"],
                      UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                      UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                      tSNE_1 = data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                      tSNE_2 = data@reductions$tsne@cell.embeddings[,"tSNE_2"])
  plotx <- cbind(plotx, data@meta.data)
  if(!is.null(groups)){
    for(i in 1:length(groups)){
      plotx[,groups[i]] <- data@meta.data[,groups[i]]
    }
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
                              strip_size = 15, legend_pos = "right", ncol = 2, point_size = 0.5){
  # set.seed(2022)
  p <- ggplot(plotx, aes(x = plotx[,feature1], y = plotx[,feature2], color = plotx[,color_by])) +
    geom_point(size = point_size) +
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
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
          strip.text = element_text(size = strip_size, face = "bold"),
          strip.background = element_blank())
  if(isfacet == T & !is.null(group_by)){
    p <- p + facet_wrap(~plotx[,group_by], ncol = ncol)
  }
  return(p)
}

own_violin <- function(plotx, x = "SAMPLE_ID", feature, plotx_title, col = NULL, title.size = 20, angle = 0, hjust=NULL,vjust = NULL){
  p <- ggplot(plotx, aes_string(x=x, y=feature, fill = x)) +
    geom_violin(trim=TRUE) + scale_fill_manual(values = col)+
    theme_classic()+
    # scale_y_continuous(trans='log10') +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = angle, size = ifelse(nchar(as.character(unique(plotx$SAMPLE_ID))) > 20, 12,15), hjust=hjust,vjust = vjust),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = title.size, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 15, margin=margin(10,0,0,0),face = "bold")) +
    xlab("SAMPLE_ID") + ylab("") + ggtitle(plotx_title)
  if(max(plotx[,feature]) > 0){
    p <- p + geom_jitter(size = 0.05)
  }

  return(p)
}

own_feature <- function(plotx, feature1, feature2, title_name, col, xlabel, ylabel){
  p <- ggplot(plotx, aes(plotx[,feature1], plotx[,feature2], color = SAMPLE_ID))+
    geom_point()+theme_classic()+
    scale_color_manual(values = gen_colors(col, length(unique(plotx[,"SAMPLE_ID"])))) +
    stat_cor() +
    xlab(xlabel)+ylab(ylabel)+ # ggtitle(title_name)+
    theme(legend.position = "none", axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 15, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 15, margin=margin(0,10,0,0)),
          plot.margin=unit(c(1,1,1,1),"cm"))
  # plot.title = element_text(size = title.size, face = "bold", hjust = 0.5))
  return(p)
}

own_2d_scatter <- function(current_data, reduction_method, split_var, plot_title){
  p <- DimPlot(current_data, reduction = reduction_method,split.by = split_var, cols = color_conditions$bright) +
    ggtitle(paste(plot_title, sep = "")) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.2, "cm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 15, margin=margin(0,10,0,0)))
  return(p)
}

own_umap_metrics <- function(current_data, metrics,plot_title, sample_name){
  plotx <- current_data@meta.data
  colnames(plotx)[grep("^orig.ident$", colnames(plotx), ignore.case = T)] <- "SAMPLE_ID"
  plotx <- plotx[,grep("old.ident", colnames(plotx), ignore.case = T, invert = T)]
  plotx$UMAP_1 <- current_data@reductions$umap@cell.embeddings[,"UMAP_1"]
  plotx$UMAP_2 <- current_data@reductions$umap@cell.embeddings[,"UMAP_2"]

  melt_plotx <- melt(plotx, id.vars = c("SAMPLE_ID", "UMAP_1", "UMAP_2"), measure.vars = metrics)
  plotx <- NULL
  temp <- NULL
  for(k in 1:length(metrics)){
    temp <- melt_plotx[which(melt_plotx$variable == metrics[k]),]
    temp <- temp[which(temp$value > quantile(temp$value, 0.1)),]
    plotx <- rbind(plotx, temp)
  }

  colnames(plotx) <- c("SAMPLE_ID", "UMAP_1", "UMAP_2", "FEATURE","COUNT")
  p <- NULL
  for(k in 1:length(metrics)){
    current <- plotx[which(plotx$FEATURE == metrics[k]),]
    if(k == 1){
      p <- ggplot(current, aes(UMAP_1, UMAP_2, color = COUNT)) +
        geom_point(size = 0.4)+
        scale_color_gradient(low = "green", high = "blue") +
        theme_classic()+ ggtitle(metrics[k])+
        theme(axis.text.x = element_text(size = 25),
              axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              # legend.key.size = unit(1, "cm"),
              plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
    }else{
      p <- p + ggplot(current, aes(UMAP_1, UMAP_2, color = COUNT)) +
        geom_point(size = 0.4)+
        scale_color_gradient(low = "green", high = "blue") +
        theme_classic()+ ggtitle(metrics[k])+
        theme(axis.text.x = element_text(size = 25),
              axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              # legend.key.size = unit(1, "cm"),
              plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
    }
  }
  p + plot_layout(ncol = 2) +
    plot_annotation(title = paste(plot_title," - ",sample_name, sep = ""),
                    theme = theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5)))
  return(p)
}

deanalysis <- function(data, current_clusters, plot_title,group=NULL,
                       de_analysis = "findallmarkers", n=10, cluster_name = "seurat_clusters"){

  DefaultAssay(data) <- "RNA"
  Idents(data) <- cluster_name
  out <- NULL
  current_data_markers <- NULL
  p_thresh <- 0.05
  # fc_threshold <- 1.25
  current_clusters <- sort(as.numeric(as.character(unique(current_clusters))), decreasing = F)
  top1 <- NULL
  topn <- NULL

  if(toupper(de_analysis) == toupper("findallmarkers")){
    de_type <- "TOP_DEGENES_IN_CLUSTERS"
    de_name <- ""
    current_data_markers <- FindAllMarkers(data, min.pct = 0.5, logfc.threshold = 0.25)
    current_data_markers <- current_data_markers[current_data_markers$p_val_adj < p_thresh,]
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]

    wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
    # temp <- current_data_markers[current_data_markers[,wt] > fc_threshold,]
    top1 <- choose_topn(current_data_markers, wt = wt, group = "cluster", n = 1, max = T)
    # wt <- colnames(current_data_markers)[grep("p_val_adj", colnames(current_data_markers), ignore.case = T)]
    # top1 <- rbind(top1,choose_topn(current_data_markers, wt = wt, group = "cluster", n = 1, max = F))
    top1 <- top1[which(!is.na(top1$gene)),]
    top1 <- unique(top1)

    wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
    topn <- choose_topn(current_data_markers, wt = wt, group = "cluster", n = n, max = T)
    # wt <- colnames(current_data_markers)[grep("p_val_adj", colnames(current_data_markers), ignore.case = T)]
    # topn <- rbind(topn,choose_topn(current_data_markers, wt = wt, group = "cluster", n = n, max = F))
    # topn <- unique(topn)

    current_data_markers <- current_data_markers[order(current_data_markers$cluster, decreasing = F),]
    top1 <- top1[order(top1$cluster, decreasing = F),]
    top1 <- unique(top1)
    topn <- topn[order(topn$cluster, decreasing = F),]
    topn <- unique(topn)

  }else if(toupper(de_analysis) == toupper("findconservedmarkers")){
    de_type <- "TOP_CONSERVED_GENES_GROUPS_IN_CLUSTERS"
    de_name <- "CONSERVED MARKERS: "

    for(i in 1:length(current_clusters)){
      if(min(table(data@meta.data[,group][Idents(data) == current_clusters[i]])) > 3){
        print(i)
        current <- FindConservedMarkers(data, ident.1 = current_clusters[i], grouping.var = group, verbose = FALSE)
        current$gene <- row.names(current)
        current$cluster <- current_clusters[i]
        current <- current[current$minimump_p_val < p_thresh,]
        if(nrow(current) > 0){
          current_data_markers <- rbind(current_data_markers, current)
          current_groups <- colnames(current)[grep("avg_log.*FC", colnames(current), ignore.case = T)]
          top1 <- rbind(top1,current[unique(apply(current[,current_groups],2,function(x){which.max(x)})),])
          for(m in 1:length(current_groups)){
            topn <- rbind(topn,current[order(current[,current_groups[m]], decreasing = T),][1:n,])
          }
          # current_groups <- colnames(current)[grep("_p_val_adj", colnames(current), ignore.case = T)]
          # top1 <- rbind(top1,current[unique(apply(current[,current_groups],2,function(x){which.min(x)})),])
          # for(m in 1:length(current_groups)){
          #   topn <- rbind(topn,current[order(current[,current_groups[m]], decreasing = F),][1:n,])
          # }
          top1 <- unique(top1)
          topn <- unique(topn)
        }
      }
    }
    current_data_markers <- current_data_markers[order(current_data_markers$minimump_p_val, decreasing = F),]
    top1 <- top1[order(top1$cluster, decreasing = F),]
    top1 <- unique(top1)
    topn <- topn[order(topn$cluster, decreasing = F),]
    topn <- unique(topn)

  }else if(toupper(de_analysis) == toupper("finddegroup")){
    de_type <- "TOP_DE_GENES_GROUPS_WITHIN_EACH_CLUSTER"
    de_name <- "DE BETWEEN GROUPS: "

    for(i in 1:length(current_clusters)){
      current <-  subset(data, idents = current_clusters[i])
      if(min(table(current$GROUP)) > 3){
        Idents(current) <- group
        current_groups <- unique(current@meta.data[,group])
        for(l in 1:length(current_groups)){
          for(m in 1:length(current_groups)){
            if(l != m){
              current_result <- FindMarkers(current, ident.1 = current_groups[l], ident.2 = current_groups[m])
              current_result <- current_result[current_result$p_val_adj < p_thresh,]
              current_result <- current_result[order(current_result$p_val_adj, decreasing = F),]
              current_result$gene <- row.names(current_result)
              current_result$cluster <- current_clusters[i]
              current_result$group1 <- current_groups[l]
              current_result$group2 <- current_groups[m]

              wt <- colnames(current_result)[grep("log.*FC", colnames(current_result), ignore.case = T)]
              top1 <- rbind(top1, current_result[which.max(current_result[,wt]),])
              # wt <- colnames(current_result)[grep("p_val_adj", colnames(current_result), ignore.case = T)]
              # top1 <- rbind(top1, current_result[which.min(current_result[,wt]),])

              wt <- colnames(current_result)[grep("log.*FC", colnames(current_result), ignore.case = T)]
              current_result <- current_result[order(current_result[,wt], decreasing = T),]
              topn <- rbind(topn, current_result[1:n,])
              # wt <- colnames(current_result)[grep("p_val_adj", colnames(current_result), ignore.case = T)]
              # current_result <- current_result[order(current_result[,wt], decreasing = F),]
              # topn <- rbind(topn, current_result[1:n,])

              top1 <- unique(top1)
              topn <- unique(topn)
              current_data_markers <- rbind(current_result, current_data_markers)
            }
          }
        }
      }
    }
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]
    top1 <- top1[order(top1$cluster, decreasing = F),]
    top1 <- unique(top1)
    topn <- topn[order(topn$cluster, decreasing = F),]
    topn <- unique(topn)
  }

  # findallmarkers
  # findconservedmarkers
  # finddegroup

  out$current_data_markers <- current_data_markers
  p <- NULL

  DefaultAssay(data) <- "RNA"
  for(k in 1:length(top1$gene)){
    if(toupper(de_analysis) %in% c(toupper("finddegroup"), toupper("findconservedmarkers"))){
      if(toupper(de_analysis) == toupper("finddegroup")){
        extra_info <- paste(top1$group1[k], " VS ", top1$group2[k], sep = "")
      }else{
        extra_info <- ""
      }
      p[[k]] <- FeaturePlot(data, features = top1$gene[k], cols = c("green", "blue"),split.by = group,
                            label = T, pt.size = 0.5, label.size = 6, max.cutoff = 'q95')+
        plot_annotation(title = paste(extra_info,"\nCLUSTER ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""),
                        theme =
                          theme(axis.text.x = element_text(size = 20),
                                axis.text.y = element_text(size = 20),
                                axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                legend.title = element_blank(),
                                legend.text = element_text(size = 15),
                                plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
      if(toupper(de_analysis) == toupper("finddegroup")){
      names(p)[k] <- paste("CLUSTER ",top1$cluster[k],": TOP GENE IN ", top1$group1[k],sep = "")
      }else{
        names(p)[k] <- top1$cluster[k]
      }
    }else{
      extra_info <- ""
      if(k == 1){
        p <- FeaturePlot(data, features = top1$gene[k], cols = c("green", "blue"),split.by = group,
                         label = T, pt.size = 0.5, label.size = 6, max.cutoff = 'q95')+
          plot_annotation(title = paste(extra_info,"\nCLUSTER ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""),
                          theme =
                            theme(axis.text.x = element_text(size = 20),
                                  axis.text.y = element_text(size = 20),
                                  axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                  axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 15),
                                  plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
      }else{
        p1 <- FeaturePlot(data, features = top1$gene[k], cols = c("green", "blue"),split.by = group,
                          label = T, pt.size = 0.5, label.size = 6, max.cutoff = 'q95')+
          plot_annotation(title = paste(extra_info,"\nCLUSTER ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""),
                          theme =
                            theme(axis.text.x = element_text(size = 20),
                                  axis.text.y = element_text(size = 20),
                                  axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                  axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 15),
                                  plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
        p <- p+p1
      }
    }
  }

  out$de_type <- de_type
  out$featureplot <- p
  out$top1 <- top1
  out$topn <- topn
  out$de_name <- de_name

  return(out)
}

choose_topn <- function(current, wt, group, n, max = T){
  current <- current[order(current[,wt], decreasing = ifelse(max == T, T, F)),]
  temp <- split(current,current[,group])
  temp <- lapply(temp, function(x){x <- x[1:n,]})
  temp <- do.call(rbind.data.frame, temp)
  return(temp)
}

group_medianexpr <- function(current_data_markers, data, ref_group = "seurat_clusters", group = "seurat_clusters", cell_type = F,n=10){

  wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
  if(cell_type == TRUE){
    current_data_markers$CELL_TYPE <- data@meta.data[match(current_data_markers$cluster, data@meta.data[,ref_group]),"CELL_TYPE"]
    top_markers <- current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = n, eval(parse(text = wt)))
  }else{
    top_markers <- current_data_markers %>% group_by(cluster) %>% top_n(n = n, eval(parse(text = wt)))
  }
  top_markers <- top_markers[order(unlist(top_markers[,grep("log.*FC", colnames(current_data_markers), ignore.case = T)]), decreasing = T),]
  top_markers$gene <- factor(top_markers$gene, levels = c(as.character(unlist(unique(top_markers[order(unlist(top_markers[,grep("log.*FC", colnames(current_data_markers), ignore.case = T)]), decreasing = T),"gene"])))))
  # top_markers[,group] <- factor(top_markers[,group], levels = sort(as.numeric(as.character(unique(top_markers[,group])), decreasing = F)))

  top_expr <- data.frame(t(GetAssayData(data)[which(row.names(GetAssayData(data)) %in% unique(top_markers$gene)),]))
  top_expr$population <- data@meta.data[match(row.names(top_expr), row.names(data@meta.data)),group]

  if(cell_type == FALSE){
    clusters <- sort(as.numeric(as.character(unique(top_expr$population))))
  }else{
    clusters <- sort(as.character(unique(top_expr$population)))
  }

  median_expr <- NULL
  cell_number <- NULL
  expr <- NULL

  for(i in 1:length(clusters)){
    current <- top_expr[which(top_expr$population == clusters[i]),grep("population",colnames(top_expr), ignore.case = T, invert = T)]
    cell_number <- c(cell_number,nrow(current))
    for (k in 1:ncol(current)){
      expr <- c(expr,median(current[,k]))
    }
    median_expr <- rbind(median_expr, expr)
    expr <- NULL
  }

  median_expr <- data.frame(t(median_expr))
  row.names(median_expr) <- colnames(top_expr)[grep("population",colnames(top_expr), ignore.case = T, invert = T)]
  colnames(median_expr) <- clusters
  median_expr <- median_expr[,colSums(median_expr) > 0]

  plot_median <- scale(median_expr)
  plot_median <- t(scale(t(plot_median)))
  out <- NULL
  out$plot_median <- plot_median
  out$top_markers <- top_markers
  out$top_expr <- top_expr
  # out$current_data_markers
  return(out)
}

group_medianexpr <- function(current_data_markers, data, ref_group = "seurat_clusters", group = "seurat_clusters", cell_type = F,n=10){

  wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
  if(cell_type == TRUE){
    current_data_markers$CELL_TYPE <- data@meta.data[match(current_data_markers$cluster, data@meta.data[,ref_group]),"CELL_TYPE"]
    top_markers <- current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = n, eval(parse(text = wt)))
  }else{
    top_markers <- current_data_markers %>% group_by(cluster) %>% top_n(n = n, eval(parse(text = wt)))
  }
  top_markers <- top_markers[order(unlist(top_markers[,grep("log.*FC", colnames(current_data_markers), ignore.case = T)]), decreasing = T),]
  top_markers$gene <- factor(top_markers$gene, levels = c(as.character(unlist(unique(top_markers[order(unlist(top_markers[,grep("log.*FC", colnames(current_data_markers), ignore.case = T)]), decreasing = T),"gene"])))))
  # top_markers[,group] <- factor(top_markers[,group], levels = sort(as.numeric(as.character(unique(top_markers[,group])), decreasing = F)))

  top_expr <- data.frame(t(GetAssayData(data)[which(row.names(GetAssayData(data)) %in% unique(top_markers$gene)),]))

  if(cell_type == FALSE){
    top_expr$population <- data@meta.data[match(row.names(top_expr), row.names(data@meta.data)),group]
    clusters <- sort(as.numeric(as.character(unique(top_expr$population))))
  }else{
    top_expr$population <- data@meta.data[match(row.names(top_expr), row.names(data@meta.data)),"CELL_TYPE"]
    clusters <- sort(as.character(unique(top_expr$population)))
  }

  median_expr <- NULL
  cell_number <- NULL
  expr <- NULL

  for(i in 1:length(clusters)){
    current <- top_expr[which(top_expr$population == clusters[i]),grep("population",colnames(top_expr), ignore.case = T, invert = T)]
    cell_number <- c(cell_number,nrow(current))
    for (k in 1:ncol(current)){
      expr <- c(expr,median(current[current[,k]>0,k]))
    }
    median_expr <- rbind(median_expr, expr)
    expr <- NULL
  }

  median_expr <- data.frame(t(median_expr))
  row.names(median_expr) <- colnames(top_expr)[grep("population",colnames(top_expr), ignore.case = T, invert = T)]
  colnames(median_expr) <- clusters
  # median_expr <- median_expr[,colSums(median_expr) > 0]

  plot_median <- scale(median_expr)
  plot_median <- t(scale(t(plot_median)))
  out <- NULL
  out$plot_median <- plot_median
  out$top_markers <- top_markers
  out$top_expr <- top_expr
  # out$current_data_markers
  return(out)
}

beforeafter_dimplot <- function(data1, data2, dim1, dim2, group, subtitle1, subtitle2, maintitle, titlesize, col, legend_position = "bottom"){
  p1 <- plot_bygroup(data1, x = dim1, y = dim2, group = group, plot_title = subtitle1,
                     col = col, annot = FALSE, legend_position = legend_position, point_size = 1)
  p2 <- plot_bygroup(data2,  x = dim1, y = dim2, group = group, plot_title = subtitle2,
                     col = col, annot = FALSE, legend_position = legend_position, point_size = 1)

  p <- p1+p2+
    plot_annotation(title = maintitle, theme = theme(plot.title = element_text(size = titlesize, face = "bold", hjust = 0.5)))

  return(p)
}

plot_pseudo <- function(data, reduction, group, label_size, plot_title, col, n_col,
                        cell_size = 1, traj_size = 1.5){

  p <- plot_cells(data,
                  reduction_method = reduction,
                  color_cells_by = group,
                  group_label_size = label_size,
                  graph_label_size = 6,
                  cell_size = cell_size,
                  # cell_stroke = I(2/2),
                  alpha = 1,
                  trajectory_graph_segment_size = traj_size,
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE)+
    scale_color_manual(values = gen_colors(col,n_col))+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 15),
          # legend.position = "right",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.2, "cm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

  return(p)
}

adjust_plot <- function(p,col,n_col,plot_title="", fill = F){
  if(fill == F){
    p <- p + scale_color_manual(values = gen_colors(col,n_col))
  }
  p <- p +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 15),
          # legend.position = "right",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.2, "cm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  return(p)
}

get_earliest_principal_node <- function(mono3data, group, group_element){
  cell_ids <- which(colData(mono3data)[,group] == group_element)
  # mono3data <- preprocess_cds(mono3data, num_dim = 100)
  closest_vertex <-
    mono3data@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(mono3data),])
  root_pr_nodes <-
    igraph::V(principal_graph(mono3data)[["UMAP"]])$name[as.numeric(names
                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

plot_3d <- function(d1, d2, d3, color_group, plot_title, n, lt,
                    current_text,t1,t2,t3,col = color_conditions$general){

  library(plotly)
  p <- plot_ly(x=d1,y=d2,z=d3,type="scatter3d", mode="markers",size = 0.2,
               colors = gen_colors(col, n),
               color=color_group,
               hoverinfo = "text",
               hovertext = current_text)%>%
    layout(scene = list(xaxis = list(title = t1),
                        yaxis = list(title = t2),
                        zaxis = list(title = t3)),
           title =paste("<b>",plot_title,"</b>",sep = ""),
           legend = list(title = list(text = lt)))
  return(p)

}

