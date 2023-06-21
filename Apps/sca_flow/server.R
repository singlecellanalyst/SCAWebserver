###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: Flow Analysis Pipeline
# Version: V1.2.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last-update date: 2023-06-21
# All Rights Reserved
###########################################################################################################################
library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=30000*1024^2)

library("Seurat")
library("shiny")
library("shinyjs")
library("shinyWidgets")
library("shinyscreenshot")
library("shinythemes")
library("shinyFiles")
library("shinydashboard")
library("shinyalert")
library("remotes")
library("cytolib")
library("flowCore")
library("robustbase")
library("SDMTools")
library("vite")
library("rhandsontable")
library("ggcyto")
library("flowStats")
library("gridExtra")
library("reshape2")
library("ggridges")
library("openCyto")
library("flowViz")
library("patchwork")
library("limma")
library("gplots")
library("plot3D")
library("Biobase")
library("DT")
library("cytofWorkflow")
library("FlowSOM")
library("ConsensusClusterPlus")
library("umap")
library("viridis")
library("flowWorkspace")
library("openCyto")
library("FLOWMAPR")
library("CytoExploreR")
library("flowDensity")
library("cowplot")
library("akmedoids")
library("gridGraphics")
library("patchwork")

source("DB/SCA_Flow_RShiny_Functions_V1.0.0.R")
color_conditions <- color_ini()

# example1 <- "DB/SCA_Flow_Example_From_10X.zip"
# example2 <- "DB/SCA_Flow_Metadata_Example.csv"

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
  
  results <- reactiveValues()
  
  observe({
    req(input$submit)
    req(input$files)
    req(values$proceed == 1)
    req(is.null(results$data))
    showModal(modalDialog("Initialising..", footer=NULL))
    inFile <- input$files
    print(input$phenodata$datapath)
    pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "FLOW", isDir = T)
    sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SAMPLE_ID)))
    names(sample_colors) <- unique(pheno_data$SAMPLE_ID)
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
    fs_data <- read.flowSet(path = input_dir, pattern = ".fcs", alter.names = TRUE, transformation = FALSE)
    sampleNames(fs_data) <- pheno_data[match(sampleNames(fs_data), pheno_data$FILE),grep("SAMPLE.*ID", colnames(pheno_data), ignore.case = T)]
    phenoData(fs_data)$name <- pheno_data[match(phenoData(fs_data)$name, pheno_data$FILE),grep("SAMPLE.*ID", colnames(pheno_data), ignore.case = T)]
    phenoData(fs_data)$GROUP <- pheno_data[match(sampleNames(fs_data), pheno_data$SAMPLE_ID),"GROUP"]
    phenoData(fs_data)$BATCH <- pheno_data[match(sampleNames(fs_data), pheno_data$SAMPLE_ID),"BATCH"]
    
    #######################################################################################################################################
    data <- NULL
    for(i in 1:length(fs_data)){
      data <- rbind(data,data.frame(SAMPLE_ID = sampleNames(fs_data)[i], fs_data[[i]]@exprs))
    }
    
    system(paste("rm -r ",cdir, sep = ""))
  
    current_names <- fs_data[[1]]@parameters@data[match(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)], 
                                                        fs_data[[1]]@parameters@data$name),"desc"]
    colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)] <- paste(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],ifelse(is.na(current_names), "", current_names), sep = "_")
    colnames(data) <- gsub("(.*)_$","\\1",colnames(data))
    data <- melt(data)
    colnames(data) <- c("SAMPLE_ID","CHANNEL","ARCSINH_COUNT")
    cofactor <- 150
    data$ARCSINH_COUNT <- asinh(data$ARCSINH_COUNT/cofactor)
    pheno_data <- pheno_data[order(pheno_data$GROUP, decreasing = F),]
    data$GROUP <- pheno_data[match(data$SAMPLE_ID,pheno_data$SAMPLE_ID),"GROUP"]
    data$SAMPLE_ID <- factor(data$SAMPLE_ID, unique(pheno_data$SAMPLE_ID))
    
    p1plots <- NULL
    channels <- as.character(unique(data$CHANNEL))
    for(j in 1:length(channels)){
      current <- data[data$CHANNEL == channels[j],]
      p1plots[[j]] <- ggplot(data[data$CHANNEL == channels[j],],
        aes(x = ARCSINH_COUNT, y = SAMPLE_ID, fill = SAMPLE_ID, color = SAMPLE_ID)) +
        theme_classic()+
        geom_density_ridges(alpha = 0.8) +
        ggtitle(paste("Channel: ", channels[j], sep = "")) +
        scale_fill_manual(values = sample_colors) +
        scale_color_manual(values = sample_colors)
      p1plots[[j]] <- adjust_theme(p1plots[[j]])
      names(p1plots)[j] <- channels[j]
    }
    removeModal()
    
    # 1: FSC-A vs SSC-A
    # 2: FSC-A vs FSC-H
    # 3: SSC-H vs SSC-W
    # 4: FSC-W VS FSC-H
    showModal(modalDialog("Running automated gating..", footer=NULL))
    gs <- NULL
    lastGate <- NULL
    p2plots <- NULL
    data_current <- NULL
    data_summary <- NULL
    k <- 1
    for(i in 1:length(fs_data)){
      cname <- sampleNames(fs_data)[i]
      print(paste("Running gating for ", cname, "..", sep = ""))
      gs[[i]] <- GatingSet(fs_data[i])
      gs[[i]] <- cyto_transform(gs[[i]], type = "arcsinh")
    if(length(colnames(fs_data[[i]])[grep("FSC.*A|SSC.*A", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnl <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
      chnlx <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("SSC.*A", colnames(fs_data), ignore.case = T)]
      local_min <- density(fs_data[[i]]@exprs[,chnl])$x[which(diff(sign(diff(density(fs_data[[i]]@exprs[,chnl])$y)))==+2)+1]
      current_ff <- gh_pop_get_data(gs[[i]])
      g <- openCyto:::.mindensity(current_ff, channels = chnl, filterId = "Exclude Debris", gate_range=c(local_min[1]-500,local_min[2]-1000))
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g)
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "root")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 0.7) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 1: Exclude Debris",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Exclude Debris"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      k <- k + 1
    }
    
    if(length(colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnl <- colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
      chnlx <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
      
      g <- openCyto:::.singletGate(fs_data[[i]], channels = chnl, filterId = "Exclude Doublets or Multiplets")
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g, parent = "Exclude Debris")
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "Exclude Debris")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 1.5) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 2: Exclude Doublets or Multiplets",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Exclude Doublets or Multiplets"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      k <- k + 1
    }
    
    if(length(colnames(fs_data[[i]])[grep("SSC.*W|SSC.*H", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnlx <- colnames(fs_data[[i]])[grep("SSC.*W", colnames(fs_data[[i]]), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("SSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
      
      g <- gate_flowclust_2d(fs_data[[i]], xChannel = chnlx, yChannel = chnly, K = 3, filterId = "Single Cells Gate")
      g@filterId <- "Single Cells Gate"
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g, parent = "Exclude Doublets or Multiplets")
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "Exclude Doublets or Multiplets")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 1.5) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 2: Exclude Doublets or Multiplets",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Single Cells Gate"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      k <- k + 1
    }
    
    if(length(colnames(fs_data[[i]])[grep("FSC.*H|FSC.*W", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnlx <- colnames(fs_data[[i]])[grep("FSC.*W", colnames(fs_data[[i]]), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
      
      g <- gate_flowclust_2d(fs_data[[i]], xChannel = chnlx, yChannel = chnly, K = 3, filterId = "Single Cells Gate 2")
      g@filterId <- "Single Cells Gate 2"
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g, parent = "Single Cells Gate")
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "Single Cells Gate")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 1.5) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 4: Single Cells Gate 2",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Single Cells Gate 2"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      k <- k + 1
    }
      
    data_current[[i]] <- gh_pop_get_data(gs[[i]][[1]], y = lastGate, inverse.transform = T)
    ps <- data.frame(gs_pop_get_count_fast(gs[[i]][[1]]))
    ps$Percent_of_parent <- ps$Count/ps$ParentCount
    data_summary <- rbind(data_summary, ps)
    }
    
    data_summary$name <- factor(data_summary$name, levels = sort(unique(data_summary$name)))
    gating_colors <- gen_colors(color_conditions$manycolors, length(unique(data_summary$Population)))
    names(gating_colors) <- unique(data_summary$Population)
    
    group_colors <- gen_colors(color_conditions$general, length(unique(pheno_data$GROUP)))
    names(group_colors) <- unique(pheno_data$GROUP)
    
    batch_colors <- gen_colors(color_conditions$warm, length(unique(pheno_data$BATCH)))
    names(batch_colors) <- unique(pheno_data$BATCH)
    
    data_current <- as(data_current,"flowSet")
    sampleNames(data_current) <- sampleNames(fs_data)
    pData(data_current)$name <- sampleNames(fs_data)
    
    data <- NULL
    for(i in 1:length(data_current)){
      data <- rbind(data,data.frame(SAMPLE_ID = pData(data_current)$name[i], asinh(exprs(data_current[[i]])/cofactor)))
    }
    
    if(length(unique(pheno_data$BATCH)) > 1){
      data[,which(colnames(data) != "SAMPLE_ID")] <- apply(data[,which(colnames(data) != "SAMPLE_ID")], 2, function(x){x <- lm(x ~ pheno_data[match(data$SAMPLE_ID,pheno_data$SAMPLE_ID),"BATCH"])$residual})
    }
    
    current_names <- fs_data[[1]]@parameters@data[match(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],fs_data[[1]]@parameters@data$name),"desc"]
    colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)] <- paste(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],ifelse(is.na(current_names), "", current_names), sep = "_")
    colnames(data) <- gsub("(.*)_$","\\1",colnames(data))
    
    melt_data <- melt(data)
    colnames(melt_data) <- c("SAMPLE_ID","CHANNEL","ARCSINH_COUNT")
    melt_data$GROUP <- pheno_data[match(melt_data$SAMPLE_ID,pheno_data$SAMPLE_ID),"GROUP"]
    melt_data$SAMPLE_ID <- factor(melt_data$SAMPLE_ID, unique(pheno_data$SAMPLE_ID))
    
    channels <- as.character(unique(melt_data$CHANNEL))
    p3plots <- NULL
    for(j in 1:length(channels)){
      p3plots[[j]] <- ggplot(melt_data[melt_data$CHANNEL == channels[j],], aes(x = ARCSINH_COUNT, y = SAMPLE_ID, fill = SAMPLE_ID, color = SAMPLE_ID)) +
        geom_density_ridges(alpha = 0.8) +
        theme_classic()+
        ggtitle(paste("Channel: ", channels[j], sep = "")) +
        scale_fill_manual(values = sample_colors) +
        scale_color_manual(values = sample_colors)
      p3plots[[j]] <- adjust_theme(p3plots[[j]])
      names(p3plots)[j] <- channels[j]
    }

    data_summary$Batch <- pheno_data[match(data_summary$name, pheno_data$SAMPLE_ID),"BATCH"]
    p4plots <- data_summary
    
    current <- split(data_summary, data_summary$name)
    current <- lapply(current, function(x){
      current <- data.frame(name = unique(x$name),Total_Cell = max(x$Count),
                            Filtered_Cell_Count = min(x$Count))
    })
    current <- do.call(rbind.data.frame, current)
    pheno_data$Total_Cell <- current[match(pheno_data$SAMPLE_ID, current$name),"Total_Cell"]
    pheno_data$Filtered_Cell_Count <- current[match(pheno_data$SAMPLE_ID, current$name),"Filtered_Cell_Count"]
    
    removeModal()
    
    results$data <- data
    results$pheno_data <- pheno_data
    results$data_summary <- data_summary
    results$p1plots <- p1plots
    results$p2plots <- p2plots
    results$p3plots <- p3plots
    results$p4plots <- p4plots
    results$sample_colors <- sample_colors
    results$gating_colors <- gating_colors
    results$group_colors <- group_colors
    results$batch_colors <- batch_colors
    results$project_name <- project_name
    #######################################################################################################################################
    channels <- sort(names(results$p1plots))
    annot_names <- sort(unique(results$data$SAMPLE_ID))
    results$channels <- channels
    results$annot_names <- annot_names
    
    markers <- colnames(data)[grep("SampleID|Sample.ID|SAMPLE_ID|FSC.*W|FSC.*A|FSC.*H|SSC.*A|SSC.*H|SSC.*W|Time|Hoechst|^[A-Z]+[0-9]+DI$|bead|L\\/D|Dead|Live",colnames(data), ignore.case = T, invert = T)]
    markers <- sort(markers)
    results$allmarkers <- markers
    
    removeModal()
    
    updateSelectInput(session, inputId = 'p1id', label = 'Choose a channel to display', choices = channels, selected = channels[1])
    updateSelectInput(session, inputId = 'p21id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p22id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p23id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p24id', label = 'Choose a sample to display', choices = annot_names, selected = annot_names[1])
    updateSelectInput(session, inputId = 'p3id', label = 'Choose a channel to display', choices = channels, selected = channels[1])
    updateSelectInput(session, inputId = 'p4id', label = 'Choose a category to display', choices = c("Population","Batch"), selected = "Population")
    
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
    
    if(length(input$selectmarkers) < 6 & input$submitmarkers){
      showModal(modalDialog("Please select at least five markers before you proceed", footer=NULL, easyClose = T))
    }
    
    if(length(input$selectmarkers) >= 5 & input$submitmarkers){
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
 
      results$data <- results$data[,grep("FSC.*W|FSC.*A|FSC.*H|SSC.*A|SSC.*H|SSC.*W|Time|Hoechst|^[A-Z]+[0-9]+DI$|bead|L\\/D|Dead|Live", colnames(results$data), ignore.case = T, invert = T)]
      print("results$data initial:")
      print(head(results$data))
      print("which(toupper(colnames(results$data)) %in% c(SAMPLE_ID,input$selectmarkers)):")
      print(which(toupper(colnames(results$data)) %in% c("SAMPLE_ID",input$selectmarkers)))
      results$data <- results$data[,which(toupper(colnames(results$data)) %in% c("SAMPLE_ID",toupper(input$selectmarkers)))]
      results$selectmarkers <- input$selectmarkers
      
      removeModal()
      
      showModal(modalDialog("Running samplewise median expression and distance..", footer=NULL))
      
      files <- unique(results$data$SAMPLE_ID)
      files
      print("files:")
      print(files)
      
      plot_median <- NULL
      expr <- NULL
      
      print("results$data:")
      print(head(results$data))
      
      for(i in 1:length(files)){
        current <- results$data[which(results$data$SAMPLE_ID == files[i]),]
        current <- current[,grep("SampleID|Sample.ID|SAMPLE_ID|Time|file|Hoechst|^[A-Z]+[0-9]+DI$|bead|L\\/D|Dead|Live",colnames(current), invert = T, ignore.case = T)]
        print("colnames(current):")
        print(colnames(current))
        print(head(current))
        for (j in 1:ncol(current)){
          current[,j] <- as.numeric(as.character(current[,j]))
          expr <- c(expr,median(current[current[,j] > 0,j]))
        }
        plot_median <- rbind(plot_median, expr)
        expr <- NULL
      }
      
      plot_median <- data.frame(t(plot_median))
      row.names(plot_median) <- colnames(results$data[,grep("SAMPLE.*ID", colnames(results$data), ignore.case = T, invert = T)])
      colnames(plot_median) <- files
      
      p5plots <- NULL
      p6plots <- NULL
      
      if(ncol(plot_median) > 2){
        mds <- plotMDS(plot_median, plot = FALSE)
        pca_out <- prcomp(t(plot_median), center = TRUE, scale. = FALSE)
        ggdf <- data.frame(SAMPLE_ID = colnames(plot_median), MDS1 = mds$x, MDS2 = mds$y, PC1 = pca_out$x[,1], PC2 = pca_out$x[,2])
        ggdf$GROUP <- results$pheno_data[match(ggdf$SAMPLE_ID,results$pheno_data$SAMPLE_ID),"GROUP"]
        ggdf$CELL_COUNT <- results$pheno_data[match(ggdf$SAMPLE_ID, results$pheno_data$SAMPLE_ID),"Filtered_Cell_Count"]
        
        p5plots <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = GROUP, size = CELL_COUNT, label = SAMPLE_ID)) +
          geom_point(alpha = 0.8) +
          scale_size_continuous(range = c(6, 12))+
          geom_label_repel(show.legend = F) +
          # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
          theme_bw() +
          ggtitle("Multidimensional Scaling Plot") +
          scale_color_manual(values = results$group_colors)
        p5plots <- adjust_theme(p5plots)
        
        p6plots <- ggplot(ggdf, aes(x = PC1, y = PC2, color = GROUP, size = CELL_COUNT, label = SAMPLE_ID)) +
          geom_point(alpha = 0.8) +
          scale_size_continuous(range = c(6, 12))+
          geom_label_repel(show.legend = F) +
          # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
          theme_bw() +
          ggtitle("Principle Component Analysis Plot") +
          scale_color_manual(values = results$group_colors)
        p6plots <- adjust_theme(p6plots)
      }
      
      p7data <- plot_median
      p7data[is.na(p7data)] <- 0
      p7data <- (scale((p7data)))
      p7data <- t(scale(t(p7data)))
      
      print("p7data:")
      print(head(p7data))
      
      p8data <- plot_median
      colnames(p8data) <- paste(colnames(p8data), results$pheno_data[match(colnames(p8data), results$pheno_data$SAMPLE_ID),"GROUP"])
      p8data <- as.dendrogram(hclust(as.dist(1-cor((p8data)))))
      
      NRS <- function(x, ncomp = 3){
        pr <- prcomp(x, center = TRUE, scale. = FALSE)
        score <- rowSums(outer(rep(1, ncol(x)),pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
        return(score)
      }
      
      nrs_sample <- NULL
      current <- results$data[,grep("SAMPLE.*ID",colnames(results$data), ignore.case = T, invert = T)]
      for(i in 1:length(unique(results$pheno_data$SAMPLE_ID))){
        nrs_sample <- rbind(nrs_sample, NRS(current[which(results$data$SAMPLE_ID == unique(results$pheno_data$SAMPLE_ID)[i]),]))
      }
      
      rownames(nrs_sample) <- unique(results$pheno_data$SAMPLE_ID)
      nrs_sample <- data.frame(nrs_sample)
      nrs <- colMeans(nrs_sample, na.rm = TRUE)
      markers_ord <- names(sort(nrs, decreasing = TRUE))
      nrs_sample$SAMPLE_ID <- rownames(nrs_sample)
      ggdf <- melt(nrs_sample, id.var = "SAMPLE_ID",
                   value.name = "nrs", variable.name = "Markers")
      colnames(ggdf) <- c("SAMPLE_ID","Markers","NRScore")
      ggdf$Markers <- factor(ggdf$Markers, levels = markers_ord)
      ggdf <- merge(ggdf, results$pheno_data, by.x = "SAMPLE_ID", by.y = "SAMPLE_ID")
      
      p9plots <- ggplot(ggdf, aes(x = Markers, y = NRScore)) + 
        # geom_point(aes(color = SAMPLE_ID), alpha = 0.8,
        # position = position_jitter(width = 0.3, height = 0)) +
        geom_boxplot(aes(fill = GROUP), alpha = 0.8, outlier.color = NA) +
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
        scale_color_manual(values = results$sample_colors)+
        scale_fill_manual(values = results$group_colors)+
        ylab("Non-Redundancy Score (NRScore)")
      p9plots <- adjust_theme(p9plots, xangle = 45, hejust = 1, vejust = 1)
      
      samples <- unique(results$data$SAMPLE_ID)
      data0 <- NULL
      
      for(i in 1:length(samples)){
        data0[[i]] <- results$data[which(results$data$SAMPLE_ID == samples[i]),grep("SAMPLE_ID|SAMPLEID|SAMPLE.ID",colnames(results$data), ignore.case = T, invert = T)]
      }
      
      for (i in 1:length(samples)){
        meta <- data.frame(name=colnames(data0[[i]]),desc=colnames(data0[[i]]))
        meta$range <- apply(apply(data0[[i]],2,range),2,diff)
        meta$minRange <- apply(data0[[i]],2,min)
        meta$maxRange <- apply(data0[[i]],2,max)
        data0[[i]] <- new("flowFrame",exprs=as.matrix(data0[[i]]),parameters=AnnotatedDataFrame(meta))
      }
      removeModal()
      
      showModal(modalDialog("Estimating an optimal cluster..", footer=NULL))
      
      fs_data = as(data0,"flowSet")
      pData(fs_data)$name <- samples
      som_input <- ReadInput(fs_data)
      set.seed(59)
      som <- BuildSOM(som_input)
      codes <- som$map$codes
      print("codes:")
      print(head(codes))
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
      ck2 <- min(PAC_turn[PAC_turn>3]) + 1
      if(ck1 != ck2){
        chosen_k <- min(ck1, ck2)
      }else{
        chosen_k <- ck1
      }
      
      row.names(results$data) <- paste("Cell",1:nrow(results$data), sep = "")
      PAC$Label <- ifelse(PAC$K == chosen_k, paste("Recommended K = ", chosen_k, sep = ""), "")
      
      p101plots <- ggplot(PAC, aes(x= K, y= PAC, label = Label)) + geom_line(colour = "grey") +
        ggtitle("Recommended K Number of Clusters Based on PAC Method")+
        theme_classic(base_size = 20)+ geom_text_repel(
          max.overlaps = Inf,force=1,
          point.padding = 0, # additional padding around each point
          min.segment.length = 0, # draw all line segments
          max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
          box.padding = 0.3, size = 10, colour = "red")+
        geom_point(colour = ifelse(PAC$K == chosen_k, "red", "grey"), size = ifelse(PAC$K == chosen_k, 5, 2))
      
      results$p5plots <- p5plots
      results$p6plots <- p6plots
      results$p7data <- p7data
      results$p8data <- p8data
      results$p9plots <- p9plots
      results$p101plots <- p101plots
      results$recommend_k <- chosen_k
      results$mc <- mc
      results$som <- som
      #######################################################################################################################################
      removeModal()
      
      updateSelectInput(session, inputId = 'p9id', label = 'Choose a category to display', choices = c("NRScore","PCA"), selected = "NRScore")

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
      req(is.null(results$cluster_colors))
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
        current <- results$data[which(results$data$population == pop_list[i]),grep("SAMPLE.*ID|population",colnames(results$data), ignore.case = T, invert = T)]
        cell_number <- c(cell_number,nrow(current))
        for (j in 1:ncol(current)){
          expr <- c(expr,median(current[current[,j] > 0,j]))
        }
        plot_median <- rbind(plot_median, expr)
        expr <- NULL
      }
      
      row.names(plot_median) <- paste("Cluster_", pop_list, ":",cell_number,sep = "")
      colnames(plot_median) <- toupper(gsub("^.*?_","",colnames(results$data)[grep("SAMPLE.*ID|population",colnames(results$data), ignore.case = T, invert = T)]))
      plot_median <- data.frame(plot_median)
      print("plot_median:")
      print(plot_median)
      plot_median[is.na(plot_median)] <- 0
      dista <- hclust(as.dist(1-cor(t(plot_median))), method = "complete")
      plot_median <- scale(plot_median)
      plot_median <- t(scale(t(plot_median)))
      plot_median <- as.matrix(plot_median)
      
      data_meta <- results$data[,grep("SAMPLE.*ID|population", colnames(results$data), ignore.case = T)]
      x <- data.frame(t(results$data[,grep("SAMPLE.*ID|population", colnames(results$data), ignore.case = T, invert = T)]))
      x <- CreateSeuratObject(counts = x)
      x$PROJECT <- results$project_name
      x@meta.data <- cbind(x@meta.data,data_meta)
      x$orig.ident <- x$SAMPLE_ID
      x@assays$RNA@data <- x@assays$RNA@counts
      x <- ScaleData(x)
      x <- RunPCA(x, features = row.names(x))
      x <- RunUMAP(x, reduction = "pca", dims = 1:ifelse(length(x@reductions$pca) < 30, length(x@reductions$pca), 30))
      
      plotx <- data.frame(UMAP_1 = x@reductions$umap@cell.embeddings[,"UMAP_1"],
                          UMAP_2 = x@reductions$umap@cell.embeddings[,"UMAP_2"],
                          PC_1 = x@reductions$pca@cell.embeddings[,"PC_1"],
                          PC_2 = x@reductions$pca@cell.embeddings[,"PC_2"],
                          CELLL_ID = row.names(x@meta.data))
      
      plotx$CLUSTER <- factor(results$data[match(plotx$CELLL_ID, row.names(results$data)), "population"], levels = c(unique(sort(as.numeric(as.character(results$data$population))))))
      plotx$SAMPLE_ID <- results$data$SAMPLE_ID
      
      cluster_colors <- gen_colors(color_conditions$tenx,length(unique(plotx$CLUSTER)))
      names(cluster_colors) <- levels(plotx$CLUSTER)
      
      p10plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER",
                               plot_title = paste("UMAP - CLUSTER: ",results$project_name, sep = ""),
                               col = cluster_colors, annot = F, legend_position = "right",
                               numeric = T, point_size = 0.2, label_size = 6, legendsize = 15)
      
      p11plots <- ggplot(plotx,  aes(x = UMAP_1, y = UMAP_2, color = CLUSTER)) +
        geom_point(size = 0.2, alpha = 0.7) + 
        theme_classic() + facet_wrap(~SAMPLE_ID, ncol = 6) +
        guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
        scale_color_manual(values = cluster_colors) +
        ggtitle(paste("UMAP - CLUSTER (BY SAMPLES): ",results$project_name, sep = ""))
      p11plots <- adjust_theme(p11plots)
      
      p12plots <- plot_bygroup(plotx, x = "PC_1", y = "PC_2", group = "CLUSTER",
                               plot_title = paste("PCA - CLUSTER: ",results$project_name, sep = ""),
                               col = cluster_colors, annot = F, legend_position = "right",
                               numeric = T, point_size = 0.2, label_size = 6, legendsize = 15)
      p13data <- plot_median
      
      marker_list <- unique(colnames(results$data[,grep("SAMPLE.*ID|population", colnames(results$data), ignore.case = T, invert = T)]))
      melt_data <- melt(data.frame(results$data,UMAP_1 = plotx$UMAP_1, UMAP_2 = plotx$UMAP_2), 
                        id.vars = c("UMAP_1","UMAP_2", "SAMPLE_ID","population"))
      colnames(melt_data) <- c("UMAP_1","UMAP_2","SAMPLE_ID","CLUSTER","Marker","Asinh_Expression")
      
      p14plots <- ggplot(melt_data,  aes(x = UMAP_1, y = UMAP_2, color = Asinh_Expression)) +
        facet_wrap(~Marker) +
        geom_point(size = 0.2) + theme_bw()+
        scale_color_gradientn(colors = gen_colors(c("blue","cyan","green","yellow","orange","red","red4"), 100))
      p14plots <- adjust_theme(p14plots)
      
      total_count <- data.frame(table(results$data$population))
      current <- data.frame(table(results$data[,c("SAMPLE_ID","population")]))
      current <- current[current$Freq != 0,]
      colnames(current) <- c("SAMPLE_ID","CLUSTER","COUNT")
      current$CLUSTER_TOTAL_COUNT <- total_count[match(current$CLUSTER, total_count$Var1),"Freq"]
      current$PROPORTION <- current$COUNT/current$CLUSTER_TOTAL_COUNT
      # current$PROPORTION <- log2(current$COUNT+1)/log2(current$CLUSTER_TOTAL_COUNT+1)
      current$Label <- paste(current$COUNT,"(", signif(current$PROPORTION, digits = 3), ")", sep = "")
      
      node_proportion <- current
      print("node_proportion:")
      print(node_proportion)
      
      p15plots <- ggplot(node_proportion, aes(CLUSTER, PROPORTION, fill = SAMPLE_ID, label = Label))+
        geom_bar(stat="identity", alpha=0.8)+
        # coord_polar()+
        scale_fill_viridis(option = "A", discrete = T)+
        ggtitle(paste("Frequency of Samples in Each Cluster: ", results$project_name, "\n(Labels in each stacked element: Cell Number(Proportion))",sep = ""))+
        theme_classic()+
        geom_text(size = 5, position = position_stack(vjust = 0.5))
      
      p15plots <- adjust_theme(p15plots)
      
      results$marker_colors <- gen_colors(color_conditions$colorful, length(unique(row.names(x))))
      names(results$marker_colors) <- unique(row.names(x))
      
      cplotx <- x@reductions$pca@feature.loadings
      cplotx <- melt(cplotx)
      colnames(cplotx) <- c("Marker","PC_Components","Loadings")
      p92plots <- ggplot(cplotx, aes(x = PC_Components, y = Loadings, color = Marker, group = Marker)) + 
        geom_point(size = 5, alpha = 0.8) +
        geom_line()+
        geom_hline(yintercept = 0, linetype="dotted")+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
        scale_color_manual(values = results$marker_colors)+
        xlab("PC Components")+ylab("PC Loadings")
      p92plots <- adjust_theme(p92plots, xangle = 45, hejust = 1, vejust = 1)
      
    results$p10plots <- p10plots
    results$p11plots <- p11plots
    results$p12plots <- p12plots
    results$p13data <- p13data
    results$p14plots <- p14plots
    results$p15plots <- p15plots
    results$p92plots <- p92plots
    results$cluster_colors <- cluster_colors
    #######################################################################################################################################
    clusters <- sort(unique(results$data$population))
    removeModal()
    
    updateSelectInput(session, inputId = 'p16id', label = 'Choose a cluster to display', choices = clusters, selected = clusters[1])

    print("Completed!")
    return(results)
  })
  
  output$p1plot <- renderPlot({ #renderPlotly
    if(input$p1id != ""){
    showModal(modalDialog("Plotting figure 1..", footer=NULL))
    p <- results[['p1plots']][[input$p1id]]
    removeModal()
    print("done with")
    return(p)
    }
  }, height = 600, width = 900)
  
  observeEvent(input$p1plot, {
    screenshot(id="p1plot", filename = paste("1SCA_FLOW_ARCSINH_DENSITY_PLOT_CHANNEL_",input$p1id, sep = ""), scale = 2)
  })
  
  # output$p2plot <- renderPlot({
  #   if(input$p2id != ""){
  #   showModal(modalDialog("Plotting figure 2..", footer=NULL))
  #   p <- results[['p2plots']]
  #   p <- p[which(toupper(gsub("(.*):.*","\\1",names(p), ignore.case = T)) == toupper(input$p2id))]
  #   print("p2plots names:")
  #   print(names(p))
  #   # print(p[[1]])
  #   removeModal()
  #   # graphics.off()
  #   print(par("mar"))
  #   # par(mar=c(1,1,1,1))
  #   if(length(p) == 4){
  #    p <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],ncol = 2, nrow = 2,scale = 0.8, labels = c(results[['project_name']],input$p2id), label_size = 18, hjust = c(-0.7,-1.3))
  #   }else if(length(p) == 3){
  #    p <- plot_grid(p[[1]],p[[2]],p[[3]],ncol = 2, nrow = 2,scale = 0.8, labels = c(results[['project_name']],input$p2id), label_size = 18, hjust = c(-0.7,-1.3))
  #   }else if(length(p) == 2){
  #    p <- plot_grid(p[[1]],p[[2]],ncol = 2, nrow = 1,scale = 0.8, labels = c(results[['project_name']],input$p2id), label_size = 18, hjust = c(-0.7,-1.3))
  #   }else if(length(p) == 1){
  #    p <- plot_grid(p[[1]],ncol = 1, nrow = 1,scale = 0.8, labels = c(results[['project_name']],input$p2id), label_size = 18, hjust = c(-0.7,-1.3))
  #   }
  #   # grDevices::dev.off()
  #   # par("mar")
  #   # par(mar=c(1,1,1,1))
  #   print(p)
  #   }
  # }, height = 1000, width = 1400)
  
  # 1: FSC-A vs SSC-A
  # 2: FSC-A vs FSC-H
  # 3: SSC-H vs SSC-W
  # 4: FSC-W VS FSC-H
  
  output$p21plot <- renderPlot({
    if(input$p21id != ""){
      showModal(modalDialog("Plotting figure 2.1..", footer=NULL))
      p <- results[['p2plots']]
      p <- p[which(toupper(names(p)) == toupper(paste(input$p21id,":Exclude Debris", sep = "")))]
      if(length(p) > 0){
        print("p2.1plots names:")
        print(names(p))
      }else{
        text <- paste("\n   No FSC-A and SSC-A pair\n",
                      "       is available in this project for gating.")
        p <- ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      }
      removeModal()
      return(p)
    }
  }, height = 600, width = 900)
  
  observeEvent(input$p21plot, {
    screenshot(id="p21plot", filename = paste("2SCA_FLOW_GATING_STRATEGY_PART1_",input$p21id, sep = ""), scale = 2)
  })
  
  output$p22plot <- renderPlot({
    if(input$p22id != ""){
      showModal(modalDialog("Plotting figure 2.2..", footer=NULL))
      p <- results[['p2plots']]
      p <- p[which(toupper(names(p)) == toupper(paste(input$p22id,":Exclude Doublets or Multiplets", sep = "")))]
      if(length(p) > 0){
        print("p2.2plots names:")
        print(names(p))
      }else{
        text <- paste("\n   No FSC-A vs FSC-H pair\n",
                      "       is available in this project for gating.")
        p <- ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      }
      removeModal()
      return(p)
    }
  }, height = 600, width = 900)
  
  observeEvent(input$p22plot, {
    screenshot(id="p22plot", filename = paste("2SCA_FLOW_GATING_STRATEGY_PART2_",input$p22id, sep = ""), scale = 2)
  })
  
  output$p23plot <- renderPlot({
    if(input$p23id != ""){
      showModal(modalDialog("Plotting figure 2.3..", footer=NULL))
      p <- results[['p2plots']]
      p <- p[which(toupper(names(p)) == toupper(paste(input$p23id,":Single Cells Gate", sep = "")))]
      if(length(p) > 0){
        print("p2.3plots names:")
        print(names(p))
      }else{
        text <- paste("\n   No SSC-H vs SSC-W pair\n",
                      "       is available in this project for gating.")
        p <- ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      }
      removeModal()
      return(p)
    }
  }, height = 600, width = 900)
  
  observeEvent(input$p23plot, {
    screenshot(id="p23plot", filename = paste("2SCA_FLOW_GATING_STRATEGY_PART3_",input$p23id, sep = ""), scale = 2)
  })
  
  output$p24plot <- renderPlot({
    if(input$p24id != ""){
      showModal(modalDialog("Plotting figure 2.4..", footer=NULL))
      p <- results[['p2plots']]
      p <- p[which(toupper(names(p)) == toupper(paste(input$p24id,":Single Cells Gate 2", sep = "")))]
      if(length(p) > 0){
        print("p2.4plots names:")
        print(names(p))
      }else{
        text <- paste("\n   No FSC-W VS FSC-H pair\n",
                      "       is available in this project for gating.")
        p <- ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      }
      removeModal()
      return(p)
    }
  }, height = 600, width = 900)
  
  observeEvent(input$p24plot, {
    screenshot(id="p24plot", filename = paste("2SCA_FLOW_GATING_STRATEGY_PART4_",input$p24id, sep = ""), scale = 2)
  })
  
  output$p2table <- renderDT(results[['data_summary']],
                              filter = "top",
                              style="bootstrap",
                              rownames = F,
                              options = list(pageLength = 10))
  
  output$p3plot <- renderPlot({
    if(input$p3id != ""){
      showModal(modalDialog("Plotting figure 3..", footer=NULL))
      p <- results[['p3plots']][[input$p3id]]
      removeModal()
      return(p)
    }
  }, height = 600, width = 900)
  
  observeEvent(input$p3plot, {
    screenshot(id="p3plot", filename = paste("3SCA_FLOW_POST_FILTERING_ARCSINH_DENSITY_PLOT_CHANNEL_",input$p3id, sep = ""), scale = 2)
  })
  
  output$p4plot <- renderPlot({
      showModal(modalDialog("Plotting figure 4..", footer=NULL))
      plotx <- NULL
      plotx <- results[['p4plots']]
      print("input$p4id:")
      print(input$p4id)
      print(head(plotx))
      clabel <- NULL
      if(input$p4id == "Batch"){
        ccols <- results$batch_colors
        clabel <- "Batch"
      }else if(input$p4id == "Population" | input$p4id == ""){
          clabel <- "Population"
        ccols <- results$gating_colors
      }
      print("Running p4 plotting..")
      p <- NULL
      p <- ggplot(plotx, aes(x = .data[["name"]], y = .data[["Count"]], fill = .data[[clabel]], group = .data[["Population"]])) +
        geom_bar(stat="identity",alpha=0.9, size=0.5, colour = "black", position = "dodge") +
        # coord_flip() +
        theme_bw() +
        theme(legend.position = "bottom") +
        ylab("Cell Count")+xlab("Samples") +
        scale_fill_manual(values = ccols)
      p <- adjust_theme(p, legend = "bottom") + guides(fill=guide_legend(ncol=1))
      print("Finished p4 plotting..")
      removeModal()
      return(p)
  }, height = 500, width = 900)
  
  observeEvent(input$p4plot, {
    screenshot(id="p4plot", filename = paste("4SCA_FLOW_CELL_COUNT_SUMMARY_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p5plot <- renderPlot({
    showModal(modalDialog("Plotting figure 5..", footer=NULL))
    if(!is.null(results[['p5plots']])){
      p <- results[['p5plots']]
    }else{
      text <- paste("\n   More than two samples\n",
                    "     are needed in order to visualize MDS.")
      p <- ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }
    removeModal()
    return(p)
  }, height = 800, width = 1200)
  
  observeEvent(input$p5plot, {
    screenshot(id="p5plot", filename = paste("5SCA_FLOW_SAMPLE_MDS_ARCSINH_MEDIAN_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p6plot <- renderPlot({
    showModal(modalDialog("Plotting figure 6..", footer=NULL))
    if(!is.null(results[['p6plots']])){
      p <- results[['p6plots']]
    }else{
      text <- paste("\n   More than two samples\n",
                    "     are needed in order to visualize PCA")
      p <- ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }
    removeModal()
    return(p)
  }, height = 800, width = 1200)
  
  observeEvent(input$p6plot, {
    screenshot(id="p6plot", filename = paste("6SCA_FLOW_SAMPLE_PCA_ARCSINH_MEDIAN_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p7plot <- renderPlot({
    showModal(modalDialog("Plotting figure 7..", footer=NULL))
    plotx <- results[['p7data']]
    print(heatmap.2(as.matrix(plotx),margin=c(20,20), trace="none",key=T, keysize=0.8,
              dendrogram = "both",key.title = "Z-Score",
              col=jet2.col(n = 100, alpha = 1),
              srtCol=45, scale="none",
              density.info="none", cexCol=1.8,cexRow=1.2))
    removeModal()
    return(p)
    
  }, height = 600, width = 1200)
  
  observeEvent(input$p7plot, {
    screenshot(id="p7plot", filename = paste("7SCA_FLOW_HEATMAP_SAMPLE_SCALED_ARCSINH_MEDIAN_EXPRESSION_",results[['project_name']], sep = ""), scale = 2)
  })
  
  output$p8plot <- renderPlot({
    showModal(modalDialog("Plotting figure 8..", footer=NULL))
    par(mar=c(3,4,1,6))
    print(plot(results[['p8data']], horiz = TRUE))
    removeModal()
    return(p)
  }, height = 1000, width = 1000)
  
  observeEvent(input$p8plot, {
    screenshot(id="p8plot", filename = paste("8SCA_FLOW_DENDROGRAM_SAMPLES_",results[['project_name']], sep = ""), scale = 2)
  })

output$p9plot <- renderPlot({
  showModal(modalDialog("Plotting figure 9..", footer=NULL))
  print("input$p9id:")
  print(input$p9id)
  if(input$p9id == "NRScore"){
    p <- results[['p9plots']]
  }else{
    p <- results[['p92plots']]
  }
  removeModal()
  return(p)
}, height = 600, width = 900)

observeEvent(input$p9plot, {
  screenshot(id="p9plot", filename = paste("9SCA_FLOW_",input$p9id,"_ARCSINH_MARKER_",results[['project_name']], sep = ""), scale = 2)
})

output$p101plot <- renderPlot({
  showModal(modalDialog("Plotting figure PAC..", footer=NULL))
  p <- results[['p101plots']]
  removeModal()
  return(p)
}, height = 700, width = 1400)

observeEvent(input$p101plot, {
  screenshot(id="p101plot", filename = paste("101SCA_FLOW_PAC_ELBOW_PLOT_CHOSEN_CLUSTER_NUMBER_",results[['project_name']], sep = ""), scale = 2)
})

output$p10plot <- renderPlot({
  showModal(modalDialog("Plotting figure 10..", footer=NULL))
  p <- results[['p10plots']]
  removeModal()
  return(p)
}, height = 800, width = 1200)

observeEvent(input$p10plot, {
  screenshot(id="p10plot", filename = paste("10SCA_FLOW_UMAP_CLUSTERS_",results[['project_name']], sep = ""), scale = 2)
})

output$p11plot <- renderPlot({
  showModal(modalDialog("Plotting figure 11..", footer=NULL))
  p <- results[['p11plots']]
  removeModal()
  return(p)
}, height = 600, width = 1200)

observeEvent(input$p11plot, {
  screenshot(id="p11plot", filename = paste("11SCA_FLOW_UMAP_CLUSTERS_BY_SAMPLES_",results[['project_name']], sep = ""), scale = 2)
})

output$p12plot <- renderPlot({
  showModal(modalDialog("Plotting figure 12..", footer=NULL))
  p <- results[['p12plots']]
  removeModal()
  return(p)
}, height = 800, width = 1200)

observeEvent(input$p12plot, {
  screenshot(id="p12plot", filename = paste("12SCA_FLOW_PCA_CLUSTERS_",results[['project_name']], sep = ""), scale = 2)
})

output$p13plot <- renderPlot({
  showModal(modalDialog("Plotting figure 13..", footer=NULL))
  plotx <- results[['p13data']]
  print(heatmap.2(as.matrix(plotx),margin=c(10,20),key.title = "Z-Score",
                  trace="none",key=T, keysize=1,
                  dendrogram = "both",
                  col=colorRampPalette(color_conditions$orangebluecont)(100),
                  srtCol=45, scale="none",
                  density.info="none", cexCol=1.5,cexRow=1.5))
  removeModal()
  return(p)
}, height = 1000, width = 800)

observeEvent(input$p13plot, {
  screenshot(id="p13plot", filename = paste("13SCA_FLOW_HEATMAP_CLUSTERS_SCALED_ARCSINH_MEDIAN_EXPRESSION_",results[['project_name']], sep = ""), scale = 2)
})

output$p14plot <- renderPlot({
  showModal(modalDialog("Plotting figure 14..", footer=NULL))
  p <- results[['p14plots']]
  removeModal()
  return(p)
}, height = 1000, width = 1400)

observeEvent(input$p14plot, {
  screenshot(id="p14plot", filename = paste("14SCA_FLOW_UMAP_MARKER_EXPRESSIONS_",results[['project_name']], sep = ""), scale = 2)
})

output$p15plot <- renderPlot({
  showModal(modalDialog("Plotting figure 15..", footer=NULL))
  p <- results[['p15plots']]
  removeModal()
  return(p)
}, height = 1200, width = 1200)

observeEvent(input$p15plot, {
  screenshot(id="p15plot", filename = paste("15SCA_FLOW_SAMPLE_PROPORTIONS_IN_CLUSTERS_",results[['project_name']], sep = ""), scale = 2)
})

output$p16plot <- renderPlot({
  if(input$p16id != ""){
  print("input p16id:")
  print(input$p16id)
  showModal(modalDialog("Plotting figure 16..", footer=NULL))
  data <- results[['data']]
  data_1 <- data[which(data$population == input$p16id),]
  nodes <- sort(unique(data$population))
  par(mfrow=c(ifelse(ceiling(length(nodes)/4) == 1, 2, ceiling(length(nodes)/4)),4))
  par(mar=c(4,2,2,1))
  dis = rep(0,length(colnames(data)[grep("SAMPLE.*ID|population|tSNE",colnames(data),ignore.case = T,invert = T)]))
  
  for(j in grep("SAMPLE.*ID|population|tSNE",colnames(data),ignore.case = T,invert = T)){
    x=data[data[,j]>0,j]
    y=data_1[data_1[,j]>0,j]
    if (length(x)>3)
    {
      dis[j-1]=median(y)/median(x)
    }
  }
  
  for(j in grep("SAMPLE.*ID|population|tSNE",colnames(data),ignore.case = T,invert = T))
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
      text(7,2.5,paste('Cluster ',input$p16id,'(',dim(data_1)[1],')=',
                       as.character(round(median(y), digits = 1)),'',sep=''),col=2)
    }
  }
  removeModal()
  }
}, height = 500, width = 900)

observeEvent(input$p16plot, {
  screenshot(id="p16plot", filename = paste("16SCA_FLOW_MARKER_DENSITY_PLOTS_CLUSTER_",input$p16id, sep = ""), scale = 2)
})

})








