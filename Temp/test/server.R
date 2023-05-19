###########################################################################################################################
# SingleCellAnalyst.org
# Pipeline: scImmune Analysis Pipeline
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last-update date: 2023-05-15
# All Rights Reserved
###########################################################################################################################
library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=30000*1024^2)

library("igraph")

Sys.setenv("DISPLAY"=":0")

# source("DB/SCA_scImmune_RShiny_Functions_V1.0.0.R")
# color_conditions <- color_ini()

vdjdb <- readRDS(url("https://www.dropbox.com/s/80oyt32bkq6jegn/scIMMUNE_VDJDB.RDS?dl=1","rb"))
hpca.se <- readRDS(url("https://www.dropbox.com/s/2xqa8a8ussh0d2n/HumanPrimaryCellAtlasData.RDS?dl=1","rb"))

ctime <- format(Sys.time(), format = "%Y%m%d%H%M%S", tz = "Europe/Stockholm")

shinyServer(function(input, output, session) {
  values <- reactiveValues(proceed = 1)
  
  observeEvent(input$submit & input$termscheck,{
    
    if(input$termscheck == FALSE & input$submit){
      showModal(modalDialog("Please agree to our terms and conditions before you proceed", footer=NULL, easyClose = F))
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
  
  data <- reactive({
    req(input$files)
    showModal(modalDialog("Processing..depending on the number of samples as well as file sizes, time of processing may varies. Please wait patiently for results to be delivered.", footer=NULL))
    # results <- readRDS("DB/RESULTS_EXAMPLE.RDS")
    print(input$phenodata$datapath)
    pheno_data <- pheno_ini(input$phenodata$datapath, pipeline = "scIMMUNE", isDir = T)
    pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, pheno_data$CELL_TYPE, sep = "_")
    sample_colors <- gen_colors(color_conditions$manycolors, length(unique(pheno_data$SID)))
    names(sample_colors) <- unique(pheno_data$SID)
    print(input$files)
    print(input$files$datapath)
    project_name <- ifelse(input$project_name == "","My_Project",gsub("\\s+|\\(|\\)|-|\\/|\\?","",input$project_name))
    print(project_name)
    integration_method <- input$integration_method
    print("integration_method:")
    print(integration_method)
    cdir <- paste(getwd(),"/",project_name,ctime, "/",sep = "")
    dir.create(cdir)
    contig_dir <- paste(cdir, "/Contig/", sep = "")
    consensus_dir <- paste(cdir, "/Consensus/", sep = "")
    rna_dir <- paste(cdir, "/RNA/", sep = "")
    dir.create(contig_dir)
    dir.create(consensus_dir)
    dir.create(rna_dir)
    
    sample_files <- unzip(input$files$datapath, exdir = cdir)
    sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
    
    contig_pheno <- pheno_data[grep("contig", pheno_data$DATA_TYPE, ignore.case = T),]
    consensus_pheno <- pheno_data[grep("consensus", pheno_data$DATA_TYPE, ignore.case = T),]
    rna_pheno <- pheno_data[grep("rna", pheno_data$DATA_TYPE, ignore.case = T),]
    removeModal()
    
    #######################################################################################################################################
    showModal(modalDialog("Reading data..", footer=NULL))
    
    contig_list <- NULL
    for(i in 1:nrow(contig_pheno)){
      contig_list[[i]] <- sample_files[grep(contig_pheno$FILE[i],sample_files, ignore.case = T)]
      system(paste("mv \"", contig_list[[i]], "\" \"",contig_dir,"/\"", sep = ""))
      contig_list[[i]] <- read.csv(paste(contig_dir,"/",contig_pheno$FILE[i],sep = ""), header = T, sep = ",")
      names(contig_list)[i] <- contig_pheno[i,"SID"]
    }
    
    contig_underscore <- as.numeric(lapply(contig_list, function(x){ x<- length(grep("_", x$barcode))}))
    contig_hyphen <- as.numeric(lapply(contig_list, function(x){ x<- length(grep("-", x$barcode))}))
    
    if(length(which(contig_underscore == 0))>1 & length(which(contig_hyphen > 0)) == length(contig_list)){
      contig_list <- contig_list
    } else{
      for (i in seq_along(contig_list)) {
        contig_list[[i]] <- stripBarcode(contig_list[[i]], column = 1, connector = "_",
                                         num_connects = max(as.numeric(lapply(contig_list, function(x){ x<- length(grep("_", x$barcode))})))+1)
      }
    }
    
    if(nrow(rna_pheno) > 0){
      rna_list <- NULL
      for(i in 1:nrow(rna_pheno)){
        rna_list[[i]] <- sample_files[grep(rna_pheno$FILE[i],sample_files, ignore.case = T)]
        system(paste("mv \"", rna_list[[i]], "\" \"",rna_dir,"/\"", sep = ""))
      }
      rna_list <- list.files(rna_dir, full.names = T)
      rna_list <- rna_list[grep("_MACOSX",rna_list, ignore.case = T, invert = T)]
      print("rna_list:")
      print(rna_list)
    }
    
    contig_monarch <- repLoad(contig_dir,.mode = "paired")
    for(i in 1:length(contig_monarch$data)){
      names(contig_monarch$data)[i] <- pheno_data[which(unlist(lapply(gsub("(.*)\\..*","\\1",pheno_data$FILE, ignore.case = T), function(x){grepl(x,names(contig_monarch$data)[i], ignore.case = T)})) == TRUE),"SID"]
    }
    contig_monarch$meta <- pheno_data[match(names(contig_monarch$data), pheno_data$SID),]
    contig_monarch$meta$Sample <- pheno_data[match(names(contig_monarch$data), pheno_data$SID),"SID"]
    
    meta_explore <- repExplore(contig_monarch$data, .method = "volume")
    p1 <- vis(meta_explore, .by = c("Sample"), .meta = contig_monarch$meta, .test = F) + xlab("Sample")
    p2 <- vis(meta_explore, .by = c("GROUP", "CELL_TYPE"), .meta = contig_monarch$meta, .test = F)
    
    p1plots <- p1 + p2
    
    exp_len_aa <- repExplore(contig_monarch$data, .method = "len", .col = "aa")
    exp_len_nt <- repExplore(contig_monarch$data, .method = "len", .col = "nt")
    exp_cnt <- repExplore(contig_monarch$data, .method = "count")
    exp_vol <- repExplore(contig_monarch$data, .method = "volume")
    
    p1 <- vis(exp_len_aa) + facet_wrap(~Sample) + ggtitle("DISTRIBUTION OF CDR3 LENGTH: AMINO ACID")
    p11 <- vis(exp_len_nt) + facet_wrap(~Sample) + ggtitle("DISTRIBUTION OF CDR3 LENGTH: NUCLEOTIDE")
    p2 <- vis(exp_cnt) + guides(color = guide_legend(ncol = ifelse(nrow(contig_monarch$meta) > 10, 2, 1)))
    p3 <- vis(exp_vol)
    
    p2plots <- p1 + p11
    p3plots <- p2 + p3
    
    imm_pr <- repClonality(contig_monarch$data, .method = "clonal.prop")
    imm_top <- repClonality(contig_monarch$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
    imm_rare <- repClonality(contig_monarch$data, .method = "rare")
    imm_hom <- repClonality(contig_monarch$data, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
    
    p4plots <- vis(imm_top) + vis(imm_top, .by = "Sample", .meta = contig_monarch$meta)
    p5plots <- vis(imm_rare) + vis(imm_rare, .by = "Sample", .meta = contig_monarch$meta)
    p6plots <- vis(imm_hom) + vis(imm_hom, .by = c("Sample","GROUP"), .meta = contig_monarch$meta)
    
    imm_ov1 <- repOverlap(contig_monarch$data, .method = "morisita", .verbose = F, .col = "nt+v+d+j")
    imm_ov2 <- repOverlap(contig_monarch$data, .method = "morisita", .verbose = F, .col = "aa")
    
    p1 <- vis(imm_ov1, .text.size = 6) + ggtitle("CLONOTYPES SHARED\nMORISITA OVERLAP INDEX (NT+VDJ)")
    p2 <- vis(imm_ov2, .text.size = 6) + ggtitle("CLONOTYPES SHARED\nMORISITA OVERLAP INDEX (AA)")
    p7plots  <- p1 + p2
    
    # Build public repertoire table using CDR3 nucleotide sequences
    pr.ntv <- pubRep(contig_monarch$data, .col = "nt+v", .quant = "count", .verbose = F)
    pr.aav <- pubRep(contig_monarch$data, .col = "aa", .quant = "count", .verbose = F)
    gene_stats <- gene_stats()
    gene_stats <- gene_stats[grep("hs", gene_stats$alias,ignore.case = T),]
    all_genes <- colnames(gene_stats)[(!grepl("alias|species", colnames(gene_stats), ignore.case = T)) &
                                        (gene_stats != 0)]
    all_gene_usage <- NULL
    for(i in 1:length(all_genes)){
      all_gene_usage <- rbind(all_gene_usage, data.frame(Gene_Type = all_genes[i], geneUsage(contig_monarch$data, .gene = paste("hs.",all_genes[i], sep = ""))))
    }
    
    n <- 20
    groups <- unique(contig_monarch$meta$CELL_TYPE)
    names <- NULL
    p9plots <- NULL
    
    for(i in 1:length(groups)){
      for(j in 1:length(all_genes)){
        current_gu <- geneUsage(contig_monarch$data[which(toupper(names(contig_monarch$data)) %in% toupper(unlist(contig_monarch$meta[contig_monarch$meta$CELL_TYPE == groups[i],"Sample"])))], .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
        current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA|None",x, ignore.case = T))) > 0})),]
        current_gu[is.na(current_gu)] <- 0
        current_gu$Names <- gsub("NA;|None;","",current_gu$Names)
        current_gu$Names <- gsub(";NA|;None","",current_gu$Names)
        columnnames <- colnames(current_gu)
        current_gu <- split(current_gu, current_gu$Names)
        current_gu <- lapply(current_gu, function(x){
          if(nrow(x) > 1){
            x <- data.frame(Names = unique(x$Names), t(colSums(x[,grep("^Names$", colnames(x), invert = T)])))
          }else{
            x <- x
          }
        })
        current_gu <- do.call(rbind.data.frame, current_gu)
        current_gu <- as_tibble(current_gu)
        class(current_gu) <- c("immunr_gene_usage","immunr_gene_usage","tbl_df","tbl","data.frame")
        if(nrow(current_gu) > n){
          current_gu <- current_gu[order(rowMeans(current_gu[,grep("Name", colnames(current_gu), ignore.case = T, invert = T)]), decreasing = T),]
          current_gu <- current_gu[1:n,]
          current_gu <- current_gu[!is.na(current_gu$Names),]
        }else{
          n <- nrow(current_gu)
        }
        p1 <- NULL
        p1 <- vis(current_gu) + ggtitle(paste("TOP",n,": GENE USAGE: ", groups[i], ", ", toupper(all_genes[j]),
                                              "\n(All results will be plotted for those with total of <",n," hits)",sep = ""))
        if(i == 1 & j == 1){
          p9plots[[1]] <- p1
          names(p9plots)[1] <- paste(groups[i], "_", toupper(all_genes[j]), sep = "")
        }else{
          p9plots[[length(p9plots)+1]] <- p1
          names(p9plots)[length(p9plots)] <- paste(groups[i], "_", toupper(all_genes[j]), sep = "")
        }
        n <- 20
      }
    }
    
    p10plots <- NULL
    for(j in 1:length(all_genes)){
      current_gu <- geneUsage(contig_monarch$data, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
      current_gu <- current_gu[!is.na(current_gu$Names),]
      current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
      imm_gu_js <- geneUsageAnalysis(current_gu, .method = "js", .verbose = F)
      imm_gu_cor <- geneUsageAnalysis(current_gu, .method = "cor", .verbose = F)
      p1 <- vis(imm_gu_js, .title = paste("Gene usage JS-divergence: ", toupper(all_genes[j]), sep = ""), .leg.title = "JS", .text.size = 4)
      p2 <- vis(imm_gu_cor, .title = paste("Gene usage correlation: ", toupper(all_genes[j]), sep = ""), .leg.title = "Cor", .text.size = 4)
      p <- p1 + p2
      p10plots[[j]] <- p
      names(p10plots)[j] <- all_genes[j]
    }
    
    p11plots <- NULL
    for(j in 1:length(all_genes)){
      current_gu <- geneUsage(contig_monarch$data, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
      current_gu <- current_gu[!is.na(current_gu$Names),]
      current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
      row.names(current_gu) <- current_gu$Names
      # colnames(current_gu) <- gsub("(.*?)_.*","\\1",colnames(current_gu))
      p <- vis(geneUsageAnalysis(current_gu, "cosine+hclust", .verbose = F))+ggtitle(paste("Optimal number of clusters\n(Clustered based on: ", toupper(all_genes[j]), ")", sep = ""))
      print(p)
      p11plots[[j]] <- p
      names(p11plots)[j] <- all_genes[j]
    }
    
    p12plots <- NULL
    
    if(length(contig_monarch$data) > 3){
      
      for(j in 1:length(all_genes)){
        current_gu <- geneUsage(contig_monarch$data, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
        current_gu <- current_gu[!is.na(current_gu$Names),]
        current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
        imm_cl_pca <- geneUsageAnalysis(current_gu, "js+pca+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
        imm_cl_mds <- geneUsageAnalysis(current_gu, "js+mds+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
        # imm_cl_tsne <- geneUsageAnalysis(current_gu, "js+tsne+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .perp = .01, .verbose = F)
        p1 <- vis(imm_cl_pca, .plot = "clust") + ggtitle(paste("PCA: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+PCA+K-MEANS)", sep = ""))
        p2 <- vis(imm_cl_mds, .plot = "clust") + ggtitle(paste("MDS: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+MDS+K-MEANS)", sep = ""))
        # p3 <- vis(imm_cl_tsne, .plot = "clust") + ggtitle(paste("tSNE: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+tSNE+K-MEANS)", sep = ""))
        p12plots[[j]] <- p1+p2
        names(p12plots)[j] <- all_genes[j]
      }
    }
    
    p13plots <- NULL
    for(i in 1:length(contig_monarch$data)){
      p1 <- vis(spectratype(contig_monarch$data[[i]], .quant = "id", .col = "nt"))+
        ggtitle(paste("SAMPLE: ",names(contig_monarch$data)[i],sep = "")) + 
        xlab('CCDR3 Nucleotide length')
      p2 <- vis(spectratype(contig_monarch$data[[i]], .quant = "count", .col = "aa+v")) + 
        xlab('CCDR3 AA length')
      p <- p1 + p2
      p13plots[[i]] <- p
      names(p13plots)[i] <- names(contig_monarch$data)[i]
    }
    
    # Diversity estimation
    diversity_estimates <- NULL
    diversity_methods <- c("chao1", "hill", "div", "gini.simp", "inv.simp", "d50","raref")
    p14plots <- NULL
    
    for(i in 1:length(diversity_methods)){
      p14plots[[i]] <- vis(repDiversity(contig_monarch$data, .method = diversity_methods[i]))
      names(p14plots)[i] <- diversity_methods[i]
    }
    
    groups <- unique(contig_monarch$meta$CELL_TYPE)
    
    n <- 10
    p15plots <- NULL
    for(i in 1:length(contig_monarch$data)){
      tc1 <- trackClonotypes(contig_monarch$data, list(i, n), .col = "nt")
      tc2 <- trackClonotypes(contig_monarch$data, list(i, n), .col = "aa")
      p1 <- vis(tc1)+ggtitle(paste(names(contig_monarch$data)[i], ": TOP ",n," CLONOTYPES (CDR3 NT)", sep = ""))+theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=1))+
        scale_fill_manual(values = gen_colors(col = color_conditions$tableau20, n = nrow(tc1)), labels = ~ stringr::str_wrap(.x, width = 20))
      p2 <- vis(tc2)+ggtitle(paste(names(contig_monarch$data)[i], ": TOP ",n," CLONOTYPES (CDR3 AA)", sep = ""))+theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=1))+
        scale_fill_manual(values = gen_colors(col = color_conditions$tableau20, n = nrow(tc2)), labels = ~ stringr::str_wrap(.x, width = 20))
      p <- p1/p2
      print(p)
      p15plots[[i]] <- p
      names(p15plots)[i] <- names(contig_monarch$data)[i]
    }
    
    result_vdjdb <- dbAnnotate(contig_monarch$data, vdjdb, "CDR3.aa", "cdr3")
    
    # kmer
    
    p16data <- repOverlap(contig_monarch$data, .col = "v+d+j")
    ov <- scale(p16data)
    ov <- t(scale(t(ov)))
    ov[is.na(ov)] <- 0
    p17plots <- complex_heatmap(ov, col = color_conditions$BlueYellowRed, legendtitle = "Z-Score", col_title = "Repetoire Overlaps (VDJ Regions)")
    
    removeModal()
    
    showModal(modalDialog("Processing..", footer=NULL, easyClose = F))
    
    contig_types <- c("BCR","TAB","TGD")
    contig_table <- NULL
    
    bcr_list <- NULL
    
    if(length(grep("B_Cell|BCell|B.*Cell|BCR", pheno_data$CELL_TYPE, ignore.case = T)) > 0){
      bcr_list <- pheno_data[grepl("B_Cell|BCell|B.*Cell|BCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),"SID"]
    }
    
    # Very slow, package problem
    bcr_pheno <- pheno_data[grepl("B_Cell|BCell|B.*Cell|BCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),]
    if(length(bcr_list) > 0){
      contig_bcr <- combineBCR(contig_list[which(names(contig_list) %in% bcr_list)], 
                               samples = bcr_pheno[match(bcr_list,bcr_pheno$SID),"SID"],
                               removeNA = TRUE)
      for(i in 1:length(contig_bcr)){
        contig_bcr[[i]]$ID <- bcr_pheno[match(bcr_list,bcr_pheno$SID),"INDIVIDUAL_ID"][i]
      }
      names(contig_bcr) <- bcr_pheno[match(bcr_list,bcr_pheno$SID),"SID"]
      contig_table$BCR <- contig_bcr
    }else{
      contig_table$BCR <- NULL
    }
    
    removeModal()
    
    showModal(modalDialog("Processing..", footer=NULL, easyClose = F))
    
    tcr_list <- NULL
    if(length(grep("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T)) > 0){
      tcr_list <- pheno_data[grepl("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),"SID"]
      
    }
    
    tab_presence <- NULL
    tgd_presence <- NULL
    contig_tcr_tab <- NULL
    contig_tcr_tgd <- NULL
    tcr_pheno <- pheno_data[grepl("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),]
    if(length(tcr_list) > 0){
      contig_tcr <- contig_list[which(names(contig_list) %in% tcr_list)]
      tab_presence <- as.character(unlist(lapply(contig_tcr, function(x){x <- ifelse(length(grep("TRA|TRB",unique(x$chain), ignore.case = T)) > 0, "TRUE")})))
      if(length(grep("TRUE", tab_presence, ignore.case = T)) > 0){
        contig_tcr_tab <- combineTCR(contig_tcr[which(tab_presence == "TRUE")],
                                     samples = tcr_pheno[match(names(contig_tcr[which(tab_presence == "TRUE")]),tcr_pheno$SID),"SID"],
                                     removeNA = TRUE,
                                     cells = "T-AB")
        for(i in 1:length(contig_tcr_tab)){
          contig_tcr_tab[[i]]$SID <- tcr_pheno[match(names(contig_tcr[which(tab_presence == "TRUE")]),tcr_pheno$SID),"SID"][i]
        }
        contig_table$TAB <- contig_tcr_tab
      }else{
        contig_table$TAB <- NULL
      }
      
      tgd_presence <- as.character(unlist(lapply(contig_tcr, function(x){x <- ifelse(length(grep("TRG|TRD",unique(x$chain), ignore.case = T)) > 1, "TRUE", "FALSE")})))
      if(length(grep("TRUE", tgd_presence, ignore.case = T)) > 0){
        contig_tcr_tgd <- combineTCR(contig_tcr[which(tgd_presence == "TRUE")],
                                     samples = pheno_data[match(names(contig_tcr[which(tgd_presence == "TRUE")]),pheno_data$SAMPLE_ID),"SID"],
                                     removeNA = TRUE,
                                     cells = "T-GD")
        for(i in 1:length(contig_tcr_tgd)){
          contig_tcr_tgd[[i]]$ID <- pheno_data[match(names(contig_tcr[which(tgd_presence == "TRUE")]),pheno_data$SID),"SID"][i]
        }
        contig_table$TGD <- contig_tcr_tgd
      }else{
        contig_table$TGD <- NULL
      }
    }
    
    p18plots <- NULL
    p19data <- NULL
    p20plots <- NULL
    p21plots <- NULL
    p22plots <- NULL
    p23plots <- NULL
    p24plots <- NULL
    p25data <- NULL
    p26plots <- NULL
    
    k <- 1
    name_ref <- data.frame(DATA_TYPE = c("BCR","TAB","TGD"),
                           DESC = c("BCR","ALPHA-BETA TCR","GAMMA-DELTA TCR"))
    
    for(i in 1:length(contig_table)){
      
      current_type <- names(contig_table)[i]
      current_name <- name_ref[match(current_type, name_ref$DATA_TYPE),"DESC"]
      
      if(length(contig_table[[i]]) > 0){
        contig_table[[i]] <- lapply(contig_table[[i]], function(x){
          x <- data.frame(x, contig_pheno[match(unique(x$sample), contig_pheno$SID),])
        })
        
        number_contig <- quantContig(contig_table[[i]], cloneCall="gene+nt", scale = T, exportTable = T)
        
        p1 <- abundanceContig(contig_table[[i]], cloneCall = "gene+nt", scale = F) + 
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right", 
                plot.title = element_text(size = 20, face = "bold")) +
          guides(fill = guide_legend(ncol = 1))
        
        p2 <- abundanceContig(contig_table[[i]], cloneCall = "gene+nt", scale = T) + 
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right", 
                plot.title = element_text(size = 20, face = "bold")) +
          guides(fill = guide_legend(ncol = 1))
        
        p18plots[[length(p18plots)+1]] <- (p1|p2)+plot_annotation(title = paste(current_name, " CLONOTYPE ABUNDANCE (VDJC GENE + CDR3 NT)", sep = ""),theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
        names(p18plots)[length(p18plots)] <- current_name
        
        current <- abundanceContig(contig_table[[i]], cloneCall = "aa", exportTable = T)
        p19data[[length(p19data)+1]] <- current[order(current$Abundance, decreasing = T),]
        names(p19data)[length(p19data)] <- current_name
        
        p1 <- lengthContig(contig_table[[i]], cloneCall="aa", chain = "both")+
          ggtitle("CDR3 AA")+
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
                axis.text.x = element_text(angle = 0, hjust=1),
                plot.title = element_text(size = 20, face = "bold")) +
          guides(fill = guide_legend(ncol = 1))+
          scale_fill_tableau()
        p1 <- adjust_theme(p1)
        
        p2 <- lengthContig(contig_table[[i]], cloneCall="nt", chain = "both")+
          ggtitle("CDR3 NT")+
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
                axis.text.x = element_text(angle = 0, hjust=1),
                plot.title = element_text(size = 20, face = "bold"))+
          guides(fill = guide_legend(ncol = 1))+
          scale_fill_tableau()
        p2 <- adjust_theme(p2)
        
        p20plots[[length(p20plots)+1]] <- (p1/p2)+plot_annotation(title = paste(current_name," CONTIG LENGTH", sep = ""),
                                                                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
        names(p20plots)[length(p20plots)] <- current_name
        
        n <- 10
        temp <- compareClonotypes(contig_table[[i]], cloneCall="aa", exportTable = T)
        # temp <- compareClonotypes_debugged(contig_table[[i]], cloneCall = "aa")
        temp <- split(temp, temp$Sample)
        temp <- lapply(temp, function(x){
          x<- x[order(x$Proportion, decreasing = T),]
          x <- x[1:n,]})
        temp <- do.call(rbind.data.frame, temp)
        p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="aa", graph = "alluvial")
        p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15)
        p21plots[[k]] <- p21plots[[k]]+
          scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
          guides(fill = guide_legend(ncol = ceiling(length(unique(temp$Clonotypes))/12))) +
          ggtitle(paste(current_name,"\nTOP CLONOTYPES (CDR3 AA)", sep = ""))
        names(p21plots)[k] <- paste("AA_",current_name, sep = "")
        k <- k+1
        
        temp <- compareClonotypes(contig_table[[i]],cloneCall="gene", exportTable = T)
        temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
        temp <- temp[order(temp$Proportion, decreasing = T),]
        temp <- split(temp, temp$Sample)
        temp <- lapply(temp, function(x){
          x<- x[order(x$Proportion, decreasing = T),]
          x <- x[1:n,]})
        temp <- do.call(rbind.data.frame, temp)
        
        p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="gene", graph = "alluvial")
        p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15)
        p21plots[[k]] <- p21plots[[k]] + 
          scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom", 
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
          guides(fill = guide_legend(ncol = ceiling(length(unique(temp$Clonotypes))/12))) +
          ggtitle(paste(current_name,"\nTOP CLONOTYPES (VDJC GENE)", sep = ""))
        names(p21plots)[k] <- paste("GENE_",current_name, sep = "")
        k <- k+1
        
        temp <- compareClonotypes(contig_table[[i]],cloneCall="nt", exportTable = T)
        # temp <- compareClonotypes_debugged(contig_table[[i]], cloneCall = "nt")
        temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
        temp <- temp[order(temp$Proportion, decreasing = T),]
        temp <- split(temp, temp$Sample)
        temp <- lapply(temp, function(x){
          x<- x[order(x$Proportion, decreasing = T),]
          x <- x[1:n,]})
        temp <- do.call(rbind.data.frame, temp)
        
        p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="nt", graph = "alluvial")
        p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15)
        p21plots[[k]] <- p21plots[[k]] + 
          scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom", 
                axis.text.x = element_text(angle = 45, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
          guides(fill = guide_legend(ncol = 1)) +
          ggtitle(paste(current_name,"\nTOP CLONOTYPES (CDR3 NT)", sep = ""))
        names(p21plots)[k] <- paste("NT_",current_name, sep = "")
        k <- k+1
        
        temp <- compareClonotypes(contig_table[[i]],cloneCall="gene+nt", exportTable = T)
        temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
        temp <- temp[order(temp$Proportion, decreasing = T),]
        temp <- split(temp, temp$Sample)
        temp <- lapply(temp, function(x){
          x<- x[order(x$Proportion, decreasing = T),]
          x <- x[1:n,]})
        temp <- do.call(rbind.data.frame, temp)
        
        p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="gene+nt", graph = "alluvial")
        p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15)
        p21plots[[k]] <- p21plots[[k]] + 
          scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom", 
                axis.text.x = element_text(angle = 45, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
          guides(fill = guide_legend(ncol = 1)) +
          ggtitle(paste(current_name,"\nTOP CLONOTYPES (VDJC GENE + CDR3 NT)", sep = ""))
        names(p21plots)[k] <- paste("GENE_NT_",current_name, sep = "")
        k <- k+1
        
        p <- NULL
        p <- vizGenes(contig_table[[i]], gene="V", scale = T)+
          ggtitle(paste(current_name,": V GENE USAGE", sep = ""))+
          theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1, size = 10),
                plot.title = element_text(size = 15))
        p22plots[[i]] <- p
        names(p22plots)[i] <- current_name
        
        p1 <- NULL
        p2 <- NULL
        p1 <- clonalHomeostasis(contig_table[[i]], cloneCall = "gene+nt") + 
          theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45, hjust=1),
                plot.title = element_text(size = 20, face = "bold"))+
          scale_x_discrete(labels= names(contig_table[[i]]))
        p2 <- clonalProportion(contig_table[[i]], cloneCall = "gene+nt") +
          theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45, hjust=1))
        p23plots[[i]] <- p1 + p2 + plot_annotation(title = paste(current_name,": VDJC GENE + CDR3 NUCLEOTIDE CLONOTYPE", sep = ""), theme =  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
        names(p23plots)[i] <- current_name
        
        p24plots[[i]] <- clonalOverlap(contig_table[[i]], cloneCall = "gene+nt", method = "morisita") +
          scale_fill_continuous_tableau(palette = "Blue-Teal") +
          theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
                plot.title = element_text(size = 15, face = "bold")) +
          ggtitle(paste(current_name, ": MORISITA INDEX FOR CLONAL OVERLAP (VDJC GENE + CDR3 NT)", sep = ""))
        p24plots[[i]] <- adjust_theme(p24plots[[i]], title_size = 15)
        names(p24plots)[i] <- current_name
        
        t4 <- clonesizeDistribution(contig_table[[i]], cloneCall = "gene+nt", method="ward.D2", exportTable = T)
        hclust <- hclust(as.dist(t4), method = "ward.D2")
        p25data[[i]] <- as.dendrogram(hclust)
        names(p25data)[i] <- current_name # JENSEN-SHANNON DISTANCE CLUSTERING (VDJC GENE + CDR3 NT)
        
        p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "aa", group = "sample")+
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right", 
                # axis.text.x = element_text(angle = 45, hjust=1),
                plot.title = element_text(size = 10, face = "bold")) +
          guides(fill = guide_legend(ncol = 1)) +
          ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (CDR3 AA)", sep = ""))
        names(p26plots)[length(p26plots)] <- paste("AA_", current_name, sep = "")
        
        p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "gene", group = "sample")+
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right", 
                # axis.text.x = element_text(angle = 45, hjust=1),
                plot.title = element_text(size = 10, face = "bold")) +
          guides(fill = guide_legend(ncol = 1)) +
          ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (VDJC GENE)", sep = ""))
        names(p26plots)[length(p26plots)] <- paste("GENE_", current_name, sep = "")
        
        p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "nt", group = "sample")+
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right", 
                # axis.text.x = element_text(angle = 45, hjust=1),
                plot.title = element_text(size = 10, face = "bold")) +
          guides(fill = guide_legend(ncol = 1)) +
          ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (CDR3 NT)", sep = ""))
        names(p26plots)[length(p26plots)] <- paste("NT_", current_name, sep = "")
        
        p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "gene+nt", group = "sample")+
          theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right", 
                # axis.text.x = element_text(angle = 45, hjust=1),
                plot.title = element_text(size = 10, face = "bold")) +
          guides(fill = guide_legend(ncol = 1)) +
          ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (VDJC GENE + CDR3 NT)", sep = ""))
        names(p26plots)[length(p26plots)] <- paste("GENE_NT_", current_name, sep = "")
        
      }
    }
    
    removeModal()
    
    showModal(modalDialog("Running scRNA-Seq..", footer=NULL, easyClose = F))
    
    if(length(rna_list) > 0){
      data_current <- NULL
      print("rna_pheno:")
      print(nrow(rna_pheno))
      print(rna_pheno)
      for(j in 1:nrow(rna_pheno)){
        current_name <- rna_pheno[j,"SID"]
        current_file <- NULL
        current_file <- rna_list[[grep(rna_pheno[j,"FILE"], rna_list, ignore.case = T)]]
        print("current_file:")
        print(current_file)
        if(length(grep("\\.RDS$", current_file, ignore.case = T)) > 0){
          data_current[[j]] <- readRDS(current_file)
        }else if(length(grep("\\.tar.gz$|\\.gzip$|\\.tar$|\\.zip$", current_file, ignore.case = T)) > 0){
          ccdir <- NULL
          ccdir <- paste(rna_dir, gsub(".*\\/RNA\\/(.*)(\\.tar|\\.zip|\\.gzip).*","\\1",current_file), "/",sep = "")
          print("ccdir after paste:")
          print(ccdir)
          ccdir <- gsub("(.*)\\.gz$","\\1",ccdir, ignore.case = T)
          check_zip(current_file,ccdir)
          ccdir <- list.files(path = ccdir, pattern = "matrix.mtx", recursive = T, full.names = T)
          print("ccdir after list:")
          print(ccdir)
          ccdir <- gsub("(.*\\/).*","\\1",ccdir, ignore.case = T)
          print("ccdir:")
          print(ccdir)
          print(j)
          data_current[[j]] <- Read10X(ccdir)
        }else{
          data_current[[j]] <- Read10X_h5(current_file)
        }
        
        if(length(data_current[[j]]) <= 5){
          temp <- data_current[[j]]$`Gene Expression`
        }else{
          temp <- data_current[[j]]
        }
        if(length(grep("Seurat", class(data_current[[j]]), ignore.case = T)) == 0){
          data_current[[j]] <- CreateSeuratObject(counts = temp[grep("ERCC_",rownames(temp), ignore.case = T, invert = T),], project = project_name, min.cells = 3, min.features = 200)
        }
        names(data_current)[j] <- current_name
        DefaultAssay(data_current[[j]]) <- "RNA"
        Idents(data_current[[j]]) <- "orig.ident"
        data_current[[j]][["percent.mt"]] <- PercentageFeatureSet(data_current[[j]], pattern = "^MT-")
        data_current[[j]] <- subset(data_current[[j]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
        # data_current[[j]] <- suppressWarnings(SCTransform(data_current[[j]], verbose = FALSE))
        data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, rna_pheno[match(current_name, rna_pheno$SID),])
        # DefaultAssay(data_current[[j]]) <- "RNA"
        data_current[[j]] <- NormalizeData(data_current[[j]])
        data_current[[j]] <- FindVariableFeatures(data_current[[j]])
        print(paste("Complete: ", rna_pheno[j,"SID"], "..", sep = ""))
      }
      
      if(length(data_current) > 1){
        data_current <- lapply(data_current, function(x) {
          x <- ScaleData(x, verbose=FALSE, features = VariableFeatures(x))
          x <- RunPCA(x, npcs = 30, verbose = FALSE, features = VariableFeatures(x))
        })
        
        if(toupper(integration_method) == "SEURAT" | is.null(integration_method)){
          red_method <- "pca"
          integrative_features <- SelectIntegrationFeatures(object.list = data_current)
          data_anchors <- FindIntegrationAnchors(object.list = data_current,
                                                 reduction = "rpca", anchor.features = integrative_features)
          data <- IntegrateData(anchorset = data_anchors)
          DefaultAssay(data) <- "integrated"
          rm(data_anchors)
        }else if(toupper(integration_method) == "HARMONY"){
          data <- merge(data_current[[1]], data_current[c(2:length(data_current))]) # , add.cell.ids = names(data_current)
          data <- NormalizeData(data, verbose = FALSE) %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
            ScaleData(verbose = FALSE) %>% 
            RunPCA(pc.genes = data@var.genes, npcs = 30, verbose = FALSE)
          data <- RunHarmony(data, group.by.vars = "BATCH")
          red_method <- "harmony"
          DefaultAssay(data) <- "RNA"
        }
        rm(data_current)
        
      }else if(length(data_current) == 1){
        data <- data_current[[1]]
        rm(data_current)
        red_method <- "pca"
        DefaultAssay(data) <- "RNA"
      }
      
      print("data:")
      print(data)
      data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
      data <- ScaleData(data, verbose = FALSE)
      data <- RunPCA(data, verbose = FALSE)
      data <- RunUMAP(data, reduction = red_method, dims = 1:ifelse(length(data@reductions[[red_method]]) < 30, length(data@reductions[[red_method]]), 30))
      # data <- RunTSNE(data, reduction = red_method, dims = 1:ifelse(length(data@reductions[[red_method]]) < 30, length(data@reductions[[red_method]]), 30), check_duplicates = FALSE)
      data <- FindNeighbors(data, reduction = red_method, dims = 1:ifelse(length(data@reductions[[red_method]]) < 30, length(data@reductions[[red_method]]), 30))
      data <- FindClusters(data, resolution = 0.2)
      cluster_colors <- gen_colors(color_conditions$tenx, length(unique(data$seurat_clusters)))
      names(cluster_colors) <- c(sort(as.numeric(as.character(unique(data$seurat_clusters)))))
      
      Idents(data) <- "seurat_clusters"
      DefaultAssay(data) <- "RNA"
      print("Running differential expresion analysis..")
      cmarkers <- NULL
      cmarkers <- FindAllMarkers(data, min.pct = 0.5, logfc.threshold = 0.25)
      cmarkers <- cmarkers[which(cmarkers$p_val < 0.05),]
      cmarkers <- cmarkers[order(cmarkers$p_val, decreasing = F),]
      print("Running SingleR..")
      clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
                         clusters =  data$seurat_clusters,
                         ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.main)
      data$CELL_TYPE <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
      data@meta.data[which(is.na(data$CELL_TYPE)),"CELL_TYPE"] <- "Unidentifiable"
      Idents(data) <- data$CELL_TYPE
      
      celltype_colors <- gen_colors(color_conditions$monet, length(unique(data$CELL_TYPE)))
      names(celltype_colors) <- sort(as.character(unique(data$CELL_TYPE)))
      
      plotx <- gen10x_plotx(data, selected = c("PCA","UMAP"), include_meta = T)
      plotx$CLUSTER <- plotx$seurat_clusters
      # p27data <- plotx
      
      group_colors <- gen_colors(color_conditions$tableau20, length(unique(data$GROUP)))
      names(group_colors) <- sort(as.character(unique(data$GROUP)))
      
      p27plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "SID", plot_title = project_name,
                               point_size = 1, numeric = F,label_size = 10,
                               col = sample_colors, annot = F, legend_position = "right")
      p27plots <- adjust_theme(p27plots)
      
      p29plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER", plot_title = project_name,
                               point_size = 1, numeric = T,label_size = 10,
                               col = cluster_colors, annot = TRUE, legend_position = "right")
      p29plots <- adjust_theme(p29plots)
      
      p30plots <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", group_by = "GROUP",
                                    color_by = "CLUSTER", isfacet = T, xlabel = "UMAP_1", ylabel = "UMAP_2",point_size = 1,
                                    title = project_name,col = cluster_colors)
      p30plots <- adjust_theme(p30plots)
      
      p31plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = project_name,
                               point_size = 1,label_size = 10,
                               col = celltype_colors, annot = TRUE, legend_position = "right")
      p31plots <- adjust_theme(p31plots)
      
      p32plots <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", group_by = "GROUP",
                                    color_by = "CELL_TYPE", isfacet = T, xlabel = "UMAP_1", ylabel = "UMAP_2",
                                    point_size = 1,
                                    title = project_name,col = celltype_colors)
      p32plots <- adjust_theme(p32plots)
      
      orig_names <- row.names(data@meta.data)
      row.names(data@meta.data) <- paste(gsub("_CONTIG$|_SCRNA$","",data$SID, ignore.case = T),
                                         row.names(data@meta.data), sep = "_")
      row.names(data@meta.data) <- gsub("_[0-9]+$","",row.names(data@meta.data), ignore.case = T)
      cmultiplexs <- rna_pheno[grep("ALL.*CELL", rna_pheno$CELL_TYPE, ignore.case = T),"SID"]
      cmultiplex_tbcr <- pheno_data[which(pheno_data$PAIR_ID %in% unique(pheno_data[which(pheno_data$SID %in% cmultiplexs),"PAIR_ID"]) &
                                            grepl("CONTIG",pheno_data$DATA_TYPE, ignore.case = T)),"SID"]
      if(length(cmultiplexs) > 0){
        for(i in 1:length(cmultiplexs)){
          row.names(data@meta.data) <- gsub(cmultiplexs[i],gsub("_ALL.*CELL.*","",cmultiplexs[i]),row.names(data@meta.data), ignore.case = T)
        }
      }
      
      current_contigs <- NULL
      for(i in 1:length(contig_table)){
        current <- contig_table[[i]]
        for(j in 1:length(current)){
          print(paste(i,",",j,sep = ""))
          if(length(which(toupper(names(current)[j]) %in% toupper(cmultiplex_tbcr))) > 0){
            current[[j]]$barcode <- gsub("_BCR|_TCR","",current[[j]]$barcode, ignore.case = T)
          }
          current[[j]]$barcode <- gsub("_CONTIG|_CONSENSUS","",current[[j]]$barcode, ignore.case = T)
          current[[j]]$barcode <- gsub("_[0-9]+$","",current[[j]]$barcode, ignore.case = T)
          data@meta.data[which(row.names(data@meta.data) %in% current[[j]]$barcode),"DATA_TYPE"] <- unique(current[[j]][,which(toupper(colnames(current[[j]])) %in% toupper("cellType"))])[1]
        }
        current_contigs <- c(current_contigs, current)
      }
      data <- combineExpression(current_contigs, data, cloneCall="aa", proportion = T)
      # table(data$cloneType)
      data@meta.data[which(!toupper(data$DATA_TYPE) %in% toupper(c("BCR","TCR","B","T","T-AB","TGD"))),"DATA_TYPE"] <- "OTHERS"
      
      current_groups <- c("Hyperexpanded (0.1 < X <= 1)", 
                          "Large (0.01 < X <= 0.1)", 
                          "Medium (0.001 < X <= 0.01)", 
                          "Small (1e-04 < X <= 0.001)", 
                          "Rare (0 < X <= 1e-04)", NA)
      Idents(data) <- data$seurat_clusters
      plotx <- gen10x_plotx(data, selected = c("PCA","UMAP"), include_meta = T)
      # p33data <- plotx
      
      clone_cols <- c(rev(c("#1D47A0","#8BC3FA","#D1FBEC","#F3B750","#EB5A35")),"grey")
      names(clone_cols) <- current_groups
      
      p33plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "cloneType", plot_title = project_name,
                               point_size = 1, numeric = F,label_size = 10,
                               col = clone_cols, annot = "cloneType", legend_position = "right")
      p33plots <- adjust_theme(p33plots)
      
      datatype_colors <- gen_colors(color_conditions$bright, length(unique(data$DATA_TYPE)))
      names(datatype_colors) <- sort(as.character(unique(data$DATA_TYPE)))
      
      p28plots <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", group_by = "DATA_TYPE",
                                    color_by = "DATA_TYPE", isfacet = T, xlabel = "UMAP_1", ylabel = "UMAP_2",point_size = 1,
                                    title = project_name,col = datatype_colors)
      p28plots <- adjust_theme(p28plots)
      
      # current <- occupiedscRepertoire(data, x.axis = "CELL_TYPE", exportTable = T, proportion = F)
      current <- data.frame(table(data@meta.data[,c("cloneType","seurat_clusters", "GROUP", "DATA_TYPE")]))
      current <- current[current$Freq > 0,]
      current <- split(current, current$seurat_clusters)
      for(k in 1:length(current)){
        temp <- current[[k]]
        temp <- split(temp, temp$GROUP)
        for(m in 1:length(temp)){
          temp2 <- temp[[m]]
          temp2 <- split(temp2, temp2$DATA_TYPE)
          temp2 <- lapply(temp2, function(x){
            x <- data.frame(x, Prop = x$Freq/sum(x$Freq))
          })
          temp2 <- do.call(rbind.data.frame, temp2)
          temp[[m]] <- temp2
        }
        temp <- do.call(rbind.data.frame, temp)
        current[[k]] <- temp
      }
      current <- do.call(rbind.data.frame, current)
      
      print("Running occupiedscRepertoire..")
      Idents(data) <- "seurat_clusters"
      p1 <- occupiedscRepertoire(data, x.axis = "seurat_clusters", proportion = F)+
        scale_fill_manual(values = clone_cols)+
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              strip.text = element_text(size = 15),
              legend.position = "none",
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.key.size = unit(1.2, "cm"),
              axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 20),
              # axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
        scale_x_continuous(breaks=seq(min(sort(unique(as.numeric(as.character(data$seurat_clusters))))),max(sort(unique(as.numeric(as.character(data$seurat_clusters))))),1))

      p2 <- occupiedscRepertoire(data, x.axis = "seurat_clusters", proportion = T)+
        scale_fill_manual(values = clone_cols)+
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              strip.text = element_text(size = 15),
              # legend.position = "right",
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.key.size = unit(1.2, "cm"),
              axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 20),
              # axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
        scale_x_continuous(breaks=seq(min(sort(unique(as.numeric(as.character(data$seurat_clusters))))),max(sort(unique(as.numeric(as.character(data$seurat_clusters))))),1))
      
      p34plots <- p1+p2
      
      p1 <- ggplot(current, aes_string(x = "seurat_clusters", y = "Freq", fill = "cloneType")) +
        geom_bar(stat = "identity") + ylab("Single Cells") + theme_classic() + 
        theme(axis.title.x = element_blank()) + geom_text(aes(label = Freq), position = position_stack(vjust = 0.5))+
        facet_wrap(~GROUP+DATA_TYPE)+
        scale_fill_manual(values = clone_cols)+
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              strip.text = element_text(size = 15),
              legend.position = "none",
              strip.background = element_blank(),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.key.size = unit(1.2, "cm"),
              axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 20),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))
      
      p2 <- ggplot(current, aes_string(x = "seurat_clusters", y = "Prop", fill = "cloneType")) +
        geom_bar(stat = "identity") + ylab("Single Cells") + theme_classic() + 
        theme(axis.title.x = element_blank()) + geom_text(aes(label = Freq), position = position_stack(vjust = 0.5))+
        facet_wrap(~GROUP+DATA_TYPE)+
        scale_fill_manual(values = clone_cols)+
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              strip.text = element_text(size = 15),
              # legend.position = "right",
              strip.background = element_blank(),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.key.size = unit(1.2, "cm"),
              axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 20),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))
      
      p35plots <- p1+p2
      
      combined_meta <- expression2List(data, split.by = "seurat_clusters")
      combined_meta <- combined_meta[unlist(lapply(combined_meta, nrow)) > 0]
      p36plots <- NULL
      p36plots$GENE_NT <- clonalDiversity(combined_meta, cloneCall = "gene+nt") +
        ggtitle("CLONAL DIVERISTY: VDJC GENE + CDR3 NUCLEOTIDE") + 
        guides(fill=guide_legend(title="CLUSTERS")) +
        theme(plot.title = element_text(size=15, face = "bold"))
      
      p36plots$GENE <- clonalDiversity(combined_meta, cloneCall = "gene") + 
        ggtitle("CLONAL DIVERISTY: VDJC GENE") +
        guides(fill=guide_legend(title="CELL TYPES")) +
        theme(plot.title = element_text(size=15, face = "bold"))
      
      p36plots$NT <- clonalDiversity(combined_meta, cloneCall = "nt") +
        ggtitle("CLONAL DIVERISTY: CDR3 NUCLEOTIDE") + 
        guides(fill=guide_legend(title="CELL TYPES")) +
        theme(plot.title = element_text(size=15, face = "bold"))
      
      p36plots$AA <- clonalDiversity(combined_meta, cloneCall = "aa") +
        ggtitle("CLONAL DIVERISTY: CDR3 AMINO ACID") + 
        guides(fill=guide_legend(title="CELL TYPES")) +
        theme(plot.title = element_text(size=15, face = "bold"))
      
      Idents(data) <- "seurat_clusters"
      p37plots <- clonalOverlay(data, reduction = "umap", 
                                freq.cutpoint = 0, bins = 10, facet = "GROUP") + 
        guides(color = "none")+ scale_color_manual(values = cluster_colors)
      p37plots <- adjust_theme(p37plots, strip_size = 25)
      
      Idents(data) <- "CELL_TYPE"
      p38plots <- clonalOverlay(data, reduction = "umap", 
                                freq.cutpoint = 0, bins = 10, facet = "GROUP") + 
        guides(color = "none")+ scale_color_manual(values = celltype_colors)
      p38plots <- adjust_theme(p38plots, strip_size = 25)
      
      p39plots <- clonalNetwork(data,
                                reduction = "umap", 
                                identity = "seurat_clusters",
                                filter.clones = NULL,
                                filter.identity = NULL,
                                cloneCall = "aa") + scale_color_manual(values = cluster_colors)
      p39plots <- adjust_theme(p39plots)
      
      p40plots <- clonalNetwork(data,
                                reduction = "umap", 
                                identity = "CELL_TYPE",
                                filter.clones = NULL,
                                filter.identity = NULL,
                                cloneCall = "aa") + scale_color_manual(values = celltype_colors)
      p40plots <- adjust_theme(p40plots)
      
    }
    
    results <- NULL
    results$contig_monarch <- contig_monarch
    results$p1plots <- p1plots # 01.0SCA_PLOT_NUMBER_UNIQUE_CLONOTYPE_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    results$p2plots <- p2plots # 01.1SCA_PLOT_DISTR_CDR3_LENGTHS_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    # results$p3plots <- p3plots # 01.2SCA_PLOT_CLONOTYPE_ABUNDANCE_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    results$p4plots <- p4plots # 01.7SCA_PLOT_TOP_CLONAL_PROPORTIONS_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    results$p5plots <- p5plots # 01.8SCA_PLOT_RARE_CLONAL_PROPORTIONS_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    results$p6plots <- p6plots # 01.9SCA_PLOT_CLONALSPACE_HOMEOSTASIS_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    results$p7plots <- p7plots # 02.0SCA_PLOT_REPERTOIRE_OVERLAPS_ALL_SAMPLES_ANALYSIS_MIX_SAMPLES
    # results$p8plots <- p8plots
    results$p9plots <- p9plots # 02.5SCA_PLOT_TOP20_GENE_USAGE_BY_CELL_TYPES_BCR_IGHD_MIX_SAMPLES
    results$p10plots <- p10plots # 02.6SCA_PLOT_GENE_USAGE_JS-DIVERGENCE_CORRELATION_ALL_SAMPLES_ANALYSIS_ighd_MIX_SAMPLE_HEATMAP
    results$p11plots <- p11plots # 02.7SCA_PLOT_HCLUSTERING_K-MEANS_ALL_SAMPLES_ANALYSIS_ighd_MIX_SAMPLES
    results$p12plots <- p12plots # 02.8SCA_PLOT_GENE_WISE_SAMPLE_DISTANCE_PCA_MDS_IGHD_MIX_SAMPLES
    results$p13plots <- p13plots # 03.0SCA_PLOT_SPECTRATYPE_CLONALTYPE_ALL_SAMPLES_ANALYSIS_S1_BCR_MIX_SAMPLES
    results$p14plots <- p14plots # 03.1SCA_PLOT_DIVERSITY_ESTIMATES_ALL_SAMPLES_ANALYSIS_CHAO1_MIX_SAMPLES
    results$p15plots <- p15plots # 03.2SCA_PLOT_ALLUVIAL_TOP_10_MOST_ABUNDANT_CLONOTYPES_S1_BCR_MIX_SAMPLES
    # kmer
    results$p16data <- p16data # 03.7SCA_PLOT_CIRCOS_DIAGRAM_REPERTOIRE_OVERLAPS_ALL_SAMPLES_MIX_SAMPLES
    results$p17plots <- p17plots # 03.8SCA_PLOT_HEATMAP_REPERTOIRE_OVERLAPS_ALL_SAMPLES_MIX_SAMPLES
    results$p18plots <- p18plots # 04.1SCA_PLOT_CLONOTYPE_ABUNDANCE_BCR_MIX_SAMPLES
    results$p19data <- p19data # Adundance table
    results$p20plots <- p20plots # 04.3SCA_PLOT_CDR3_CONTIG_LENGTH_AA_NT_BCR_MIX_SAMPLES
    results$p21plots <- p21plots # 04.4SCA_PLOT_ALLUVIAL_CLONOTYPES_PROPORTIONS_CDR3_AA_BCR_MIX_SAMPLES
    results$p22plots <- p22plots # 04.5SCA_PLOT_CONTIG_V_GENE_USAGE_BCR_MIX_SAMPLES
    results$p23plots <- p23plots # 04.6SCA_PLOT_CLONAL_HOMEOSTASIS_PROPORTION_BCR_MIX_SAMPLES
    results$p24plots <- p24plots # 04.7SCA_PLOT_MORISITA_INDEX_CLONAL_OVERLAP_BCR_MIX_SAMPLES
    results$p25data <- p25data # 04.8SCA_PLOT_HIERARCHICAL_CLUSTERING_JENSEN_SHANNON_DISTANCE_BCR_MIX_SAMPLES
    results$p26plots <- p26plots # 04.9SCA_PLOT_CLONAL_SAMPLE_DIVERSITY_BCR_MIX_SAMPLES
    results$p27plots <- p27plots # 05.5SCA_PLOT_UMAP_BY_SAMPLE_SEURAT_V3_MIX_SAMPLES
    results$p28plots <- p28plots # 05.7SCA_PLOT_UMAP_TSNE_BY_DATA_TYPE_SEURAT_V3_MIX_SAMPLES
    results$p29plots <- p29plots # UMAP CLUSTER
    results$p30plots <- p30plots # UMAP_GROUP_CLUSTER
    results$p31plots <- p31plots # UMAP_CELL_TYPE
    results$p32plots <- p32plots # UMAP_GROUP_CELL_TYPE
    results$p33plots <- p33plots # 05.8SCA_PLOT_UMAP_BY_CLONOTYPE_SEURAT_V3_MIX_SAMPLES
    results$p34plots <- p34plots # 06.1SCA_TABLE_CLONOTYPES_BY_CELL_TYPE_CLONECALL_MIX_SAMPLES
    results$p35plots <- p35plots # 06.2SCA_TABLE_CLONOTYPES_BY_CELL_TYPE_GROUP_DATA_TYPE_CLONECALL_MIX_SAMPLES
    results$p36plots <- p36plots # 06.3SCA_PLOT_CLONO_DIVERISITY_BY_CELL_TYPE_CLONECALL_MIX_SAMPLES
    results$p37plots <- p37plots # clonalOverlay contours on clusters
    results$p38plots <- p38plots # clonalOverlay contours on celltypes
    results$p39plots <- p39plots # clonalNetwork on clusters
    results$p40plots <- p40plots # clonalNetwork on celltypes
    results$cmarkers <- cmarkers # degenes of clusters
    results$sample_colors <- sample_colors
    results$group_colors <- group_colors
    results$cluster_colors <- cluster_colors
    results$celltype_colors <- celltype_colors
    results$datatype_colors <- datatype_colors
    results$project_name <- project_name
    #######################################################################################################################################
    # saveRDS(results,"results.RDS")
    removeModal()
    
    updateSelectInput(session, inputId = 'p9id', label = 'Choose a gene to display', choices = names(p9plots), selected = names(p9plots)[1])
    updateSelectInput(session, inputId = 'p10id', label = 'Choose a gene to display', choices = names(p10plots), selected = names(p10plots)[1])
    updateSelectInput(session, inputId = 'p11id', label = 'Choose a gene to display', choices = names(p11plots), selected = names(p11plots)[1])
    updateSelectInput(session, inputId = 'p12id', label = 'Choose a gene to display', choices = names(p12plots), selected = names(p12plots)[1])
    updateSelectInput(session, inputId = 'p13id', label = 'Choose a sample to display', choices = names(p13plots), selected = names(p13plots)[1])
    updateSelectInput(session, inputId = 'p14id', label = 'Choose a diversity method', choices = names(p14plots), selected = names(p14plots)[1])
    updateSelectInput(session, inputId = 'p15id', label = 'Choose a sample to display', choices = names(p15plots), selected = names(p15plots)[1])
    updateSelectInput(session, inputId = 'p18id', label = 'Choose a cell type to display', choices = names(p18plots), selected = names(p18plots)[1])
    updateSelectInput(session, inputId = 'p19id', label = 'Choose a sample to display', choices = names(p19data), selected = names(p19data)[1])
    updateSelectInput(session, inputId = 'p20id', label = 'Choose a cell type to display', choices = names(p20plots), selected = names(p20plots)[1])
    updateSelectInput(session, inputId = 'p21id', label = 'Choose a cell type to display', choices = names(p21plots), selected = names(p21plots)[1])
    updateSelectInput(session, inputId = 'p22id', label = 'Choose a cell type to display', choices = names(p22plots), selected = names(p22plots)[1])
    updateSelectInput(session, inputId = 'p23id', label = 'Choose a cell type to display', choices = names(p23plots), selected = names(p23plots)[1])
    updateSelectInput(session, inputId = 'p24id', label = 'Choose a cell type to display', choices = names(p24plots), selected = names(p24plots)[1])
    updateSelectInput(session, inputId = 'p25id', label = 'Choose a cell type to display', choices = names(p25data), selected = names(p25data)[1])
    updateSelectInput(session, inputId = 'p26id', label = 'Choose a cell type to display', choices = names(p26plots), selected = names(p26plots)[1])
    updateSelectInput(session, inputId = 'p36id', label = 'Choose a type to display', choices = names(p36plots), selected = names(p36plots)[1])
    updateSliderInput(session, inputId = 'kmerid', label = 'Choose a k to show its top k-mers', value = 2, min = 2, max = 15, step = 1)
    updateSliderInput(session, inputId = 'kmer2id', label = 'Choose a k to show its top k-mers', value = 2, min = 2, max = 15, step = 1)
    updateSelectInput(session, inputId = 'sampleid', label = 'Choose a sample to display', choices = names(contig_monarch$data), selected = names(contig_monarch$data)[1])
    
    print("Completed!")
    return(results)
  })
  
  output$p1plot <- renderPlot({ #renderPlotly
    showModal(modalDialog("Plotting figure 1..", footer=NULL))
    p <- data()[['p1plots']]
    return(p)
    removeModal()
  }, height = 700, width = 1400)
  
  observeEvent(input$p1plot, {
    screenshot(id="p1plot", filename = paste("1SCA_scImmune_NUMBER_UNIQUE_CLONOTYPES_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p2plot <- renderPlot({
    showModal(modalDialog("Plotting figure 2..", footer=NULL))
    p <- data()[['p2plots']]
    return(p)
    removeModal()
  }, height = 700, width = 1400)
  
  observeEvent(input$p2plot, {
    screenshot(id="p2plot", filename = paste("2SCA_scImmune_DISTR_CDR3_LENGTHS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  # output$p3plot <- renderPlot({
  #   showModal(modalDialog("Plotting figure 3..", footer=NULL))
  #   p <- data()[['p3plots']]
  #   removeModal()
  #   return(p)
  # }, height = 800, width = 1400)
  # 
  # observeEvent(input$p3plot, {
  #   screenshot(id="p3plot", filename = paste("3SCA_scImmune_MDS_SAMPLE_DISTANCE_",data()[['project_name']], sep = ""), scale = 2)
  # })
  
  output$p4plot <- renderPlot({
    showModal(modalDialog("Plotting figure 4..", footer=NULL))
    p <- data()[['p4plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p4plot, {
    screenshot(id="p4plot", filename = paste("4SCA_scImmune_TOP_CLONAL_PROPORTIONS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p5plot <- renderPlot({
    showModal(modalDialog("Plotting figure 5..", footer=NULL))
    p <- data()[['p5plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p5plot, {
    screenshot(id="p5plot", filename = paste("5SCA_scImmune_RARE_CLONAL_PROPORTIONS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p6plot <- renderPlot({
    showModal(modalDialog("Plotting figure 6..", footer=NULL))
    p <- data()[['p6plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p6plot, {
    screenshot(id="p6plot", filename = paste("6SCA_scImmune_CLONALSPACE_HOMEOSTASIS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p7plot <- renderPlot({
    showModal(modalDialog("Plotting figure 7..", footer=NULL))
    p <- data()[['p7plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p7plot, {
    screenshot(id="p7plot", filename = paste("7SCA_scImmune_REPERTOIRE_OVERLAPS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  # output$p8plot <- renderPlot({
  #   showModal(modalDialog("Plotting figure 8..", footer=NULL))
  #   p <- data()[['p8plots']]
  #   removeModal()
  #   return(p)
  # }, height = 700, width = 1400)
  # 
  # observeEvent(input$p8plot, {
  #   screenshot(id="p8plot", filename = paste("8SCA_scImmune_NRS_RANKING_ALL_MARKERS_",data()[['project_name']], sep = ""), scale = 2)
  # })
  
  output$p9plot <- renderPlot({
    showModal(modalDialog("Plotting figure 9..", footer=NULL))
    p <- data()[['p9plots']][[input$p9id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p9plot, {
    screenshot(id="p9plot", filename = paste("9SCA_scImmune_GENE_USAGE_BY_CELL_TYPES_",input$p9id, sep = ""), scale = 2)
  })
  
  output$p10plot <- renderPlot({
    showModal(modalDialog("Plotting figure 10..", footer=NULL))
    p <- data()[['p10plots']][[input$p10id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p10plot, {
    screenshot(id="p10plot", filename = paste("10SCA_scImmune_GENE_USAGE_JS-DIVERGENCE_CORRELATION_",input$p10id, sep = ""), scale = 2)
  })
  
  output$p11plot <- renderPlot({
    showModal(modalDialog("Plotting figure 11..", footer=NULL))
    p <- data()[['p11plots']][[input$p11id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p11plot, {
    screenshot(id="p11plot", filename = paste("11SCA_scImmune_HCLUSTERING_K-MEANS_",input$p11id, sep = ""), scale = 2)
  })
  
  output$p12plot <- renderPlot({
    showModal(modalDialog("Plotting figure 12..", footer=NULL))
    p <- data()[['p12plots']][[input$p12id]]
    removeModal()
    return(p)
  }, height = 400, width = 900)
  
  observeEvent(input$p12plot, {
    screenshot(id="p12plot", filename = paste("12SCA_scImmune_GENE_WISE_SAMPLE_DISTANCE_PCA_MDS_",input$p12id, sep = ""), scale = 2)
  })
  
  output$p13plot <- renderPlot({
    showModal(modalDialog("Plotting figure 13..", footer=NULL))
    p <- data()[['p13plots']][[input$p13id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p13plot, {
    screenshot(id="p13plot", filename = paste("13SCA_scImmune_SPECTRATYPE_CLONALTYPE_",input$p13id, sep = ""), scale = 2)
  })
  
  output$p14plot <- renderPlot({
    showModal(modalDialog("Plotting figure 14..", footer=NULL))
    p <- data()[['p14plots']][[input$p14id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p14plot, {
    screenshot(id="p14plot", filename = paste("14SCA_scImmune_DIVERSITY_ESTIMATES_",input$p14id, sep = ""), scale = 2)
  })
  
  output$p15plot <- renderPlot({
    showModal(modalDialog("Plotting figure 15..", footer=NULL))
    p <- data()[['p15plots']][[input$p15id]]
    removeModal()
    return(p)
  }, height = 1600, width = 900)
  
  observeEvent(input$p15plot, {
    screenshot(id="p15plot", filename = paste("15SCA_scImmune_ALLUVIAL_TOP_10_MOST_ABUNDANT_CLONOTYPES_",input$p15id, sep = ""), scale = 2)
  })
  
  output$kmerplot <- renderPlot({
    showModal(modalDialog("Plotting figure 15..", footer=NULL))
    plotx <- data()[['contig_monarch']][['data']]
    kmer_length <- input$kmerid
    n <- 10
    kmers <- getKmers(plotx, kmer_length)
    kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
    kmers <- kmers[order(unlist(kmers[,2]), decreasing = T),]
    p1 <- vis(kmers, .position = "stack", .head = n)
    p2 <- vis(kmers, .head = n, .position = "dodge")
    p <- p1/p2
    removeModal()
    return(p)
  }, height = 700, width = 900)
  
  observeEvent(input$kmerplot, {
    screenshot(id="kmerplot", filename = paste("15SCA_scImmune_TOP_10_AA_KMER_SIZE_",input$kmerid, sep = ""), scale = 2)
  })
  
  output$kmer2plot <- renderPlot({
    showModal(modalDialog("Plotting figure 15..", footer=NULL))
    plotx <- data()[['contig_monarch']][['data']][[input$sampleid]]
    kmer_length <- input$kmer2id
    n <- 10
    kmers <- getKmers(plotx, kmer_length)
    kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
    kmers <- kmers[order(unlist(kmers[,2]), decreasing = T),]
    kmers_aa_stats <- kmer_profile(kmers)
    colnames(kmers_aa_stats) <- paste("KMER_POS_", 1:ncol(kmers_aa_stats), sep = "")
    p1 <- vis(kmers_aa_stats) + ggtitle(paste("POSITION FREQUENCY MATRIX: ", input$sampleid, sep = ""))+
      scale_color_manual(values = gen_colors(color_conditions$tableau20,nrow(kmers_aa_stats)))+
      theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1,vjust = 1))
    p2 <- vis(kmers_aa_stats, .plot = "seq")+scale_fill_manual(values = gen_colors(color_conditions$tenx,nrow(kmers_aa_stats)))
    p <- p1+p2
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$kmer2plot, {
    screenshot(id="kmer2plot", filename = paste("15SCA_scImmune_AA_SEQUENCE_MOTIF_KMER_SIZE_",input$kmerid, sep = ""), scale = 2)
  })
  
  output$p16plot <- renderPlot({
    showModal(modalDialog("Plotting figure 16..", footer=NULL))
    plotx <- data()[['p16data']]
    removeModal()
    vis(plotx, .plot = "circos", annotationTrack = c("grid", "axis"), preAllocateTracks = 1, grid.col = data()[['sample_colors']], transparency = 0.2)
    title(paste(data()[['project_name']],": Repertoire Overlaps (CDR3 AA)", sep = ""), cex = 15)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      circos.text(CELL_META$xcenter, 0.8, CELL_META$sector.index,cex = 2,
                  facing = "bending.outside", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA)
  }, height = 1200, width = 1200)
  
  observeEvent(input$p16plot, {
    screenshot(id="p16plot", filename = paste("16SCA_scImmune_CIRCOS_DIAGRAM_REPERTOIRE_OVERLAPS_",input$kmerid, sep = ""), scale = 2)
  })
  
  output$p17plot <- renderPlot({
    showModal(modalDialog("Plotting figure 17..", footer=NULL))
    p <- data()[['p17plots']]
    removeModal()
    return(p)
  }, height = 900, width = 1200)
  
  observeEvent(input$p17plot, {
    screenshot(id="p17plot", filename = paste("17SCA_scImmune_HEATMAP_REPERTOIRE_OVERLAPS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p18plot <- renderPlot({
    showModal(modalDialog("Plotting figure 18..", footer=NULL))
    p <- data()[['p18plots']][[input$p18id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p18plot, {
    screenshot(id="p18plot", filename = paste("18SCA_scImmune_CLONOTYPE_ABUNDANCE_",input$p18id, sep = ""), scale = 2)
  })
  
  output$p19table <- renderDT(data()[['p19data']][[input$p19id]],
                              filter = "top",
                              style="bootstrap",
                              rownames = F,
                              options = list(pageLength = 10))
  
  output$demarkers <- renderDT(data()[['cmarkers']],
                              filter = "top",
                              style="bootstrap",
                              rownames = F,
                              options = list(pageLength = 10))
  
  output$p20plot <- renderPlot({
    showModal(modalDialog("Plotting figure 20..", footer=NULL))
    p <- data()[['p20plots']][[input$p20id]]
    removeModal()
    return(p)
  }, height = 600, width = 900)
  
  observeEvent(input$p20plot, {
    screenshot(id="p20plot", filename = paste("20SCA_scImmune_CDR3_CONTIG_LENGTH_AA_NT_",input$p20id, sep = ""), scale = 2)
  })
  
  output$p21plot <- renderPlot({
    showModal(modalDialog("Plotting figure 21..", footer=NULL))
    p <- data()[['p21plots']][[input$p21id]]
    removeModal()
    return(p)
  }, height = 800, width = 900)
  
  observeEvent(input$p21plot, {
    screenshot(id="p21plot", filename = paste("21SCA_scImmune_ALLUVIAL_CLONOTYPES_PROPORTIONS_CDR3_AA_",input$p21id, sep = ""), scale = 2)
  })
  
  output$p22plot <- renderPlot({
    showModal(modalDialog("Plotting figure 22..", footer=NULL))
    p <- data()[['p22plots']][[input$p22id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p22plot, {
    screenshot(id="p22plot", filename = paste("22SCA_scImmune_CONTIG_V_GENE_USAGE_",input$p22id, sep = ""), scale = 2)
  })
  
  output$p23plot <- renderPlot({
    showModal(modalDialog("Plotting figure 23..", footer=NULL))
    p <- data()[['p23plots']][[input$p23id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p23plot, {
    screenshot(id="p23plot", filename = paste("23SCA_scImmune_CLONAL_HOMEOSTASIS_PROPORTION_",input$p23id, sep = ""), scale = 2)
  })
  
  output$p24plot <- renderPlot({
    showModal(modalDialog("Plotting figure 24..", footer=NULL))
    p <- data()[['p24plots']][[input$p24id]]
    removeModal()
    return(p)
  }, height = 600, width = 900)
  
  observeEvent(input$p24plot, {
    screenshot(id="p24plot", filename = paste("24SCA_scImmune_MORISITA_INDEX_CLONAL_OVERLAP_",input$p24id, sep = ""), scale = 2)
  })
  
  output$p25plot <- renderPlot({
    par(mar=c(3,4,1,6))
    print(plot(data()[['p25data']][[input$p25id]], horiz = TRUE))
  }, height = 500, width = 900)
  
  observeEvent(input$p25plot, {
    screenshot(id="p25plot", filename = paste("25SCA_scImmun_HIERARCHICAL_CLUSTERING_JENSEN_SHANNON_DISTANCE_",input$p25id, sep = ""), scale = 2)
  })
  
  output$p26plot <- renderPlot({
    showModal(modalDialog("Plotting figure 26..", footer=NULL))
    p <- data()[['p26plots']][[input$p26id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p26plot, {
    screenshot(id="p26plot", filename = paste("26SCA_scImmune_CLONAL_SAMPLE_DIVERSITY_",input$p26id, sep = ""), scale = 2)
  })
  
  output$p27plot <- renderPlot({
    showModal(modalDialog("Plotting figure 27..", footer=NULL))
    p <- data()[['p27plots']]
    removeModal()
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p27plot, {
    screenshot(id="p27plot", filename = paste("27SCA_scImmune_UMAP_BY_SAMPLE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p28plot <- renderPlot({
    showModal(modalDialog("Plotting figure 28..", footer=NULL))
    p <- data()[['p28plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p28plot, {
    screenshot(id="p28plot", filename = paste("28SCA_scImmune_UMAP_BY_DATA_TYPE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p29plot <- renderPlot({
    showModal(modalDialog("Plotting figure 29..", footer=NULL))
    p <- data()[['p29plots']]
    removeModal()
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p29plot, {
    screenshot(id="p29plot", filename = paste("29SCA_scImmune_UMAP_BY_CLUSTER_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p30plot <- renderPlot({
    showModal(modalDialog("Plotting figure 30..", footer=NULL))
    p <- data()[['p30plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p30plot, {
    screenshot(id="p30plot", filename = paste("30SCA_scImmune_UMAP_BY_GROUP_CLUSTER_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p31plot <- renderPlot({
    showModal(modalDialog("Plotting figure 31..", footer=NULL))
    p <- data()[['p31plots']]
    removeModal()
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p31plot, {
    screenshot(id="p31plot", filename = paste("31SCA_scImmune_UMAP_BY_CELL_TYPE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p32plot <- renderPlot({
    showModal(modalDialog("Plotting figure 32..", footer=NULL))
    p <- data()[['p32plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p32plot, {
    screenshot(id="p32plot", filename = paste("32SCA_scImmune_UMAP_BY_GROUP_CELL_TYPE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p33plot <- renderPlot({
    showModal(modalDialog("Plotting figure 33..", footer=NULL))
    p <- data()[['p33plots']]
    removeModal()
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p33plot, {
    screenshot(id="p33plot", filename = paste("33SCA_scImmune_UMAP_BY_CLONOTYPE_",input$integration_method,"_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p34plot <- renderPlot({
    showModal(modalDialog("Plotting figure 34..", footer=NULL))
    p <- data()[['p34plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p34plot, {
    screenshot(id="p34plot", filename = paste("34SCA_scImmune_CLONOTYPES_BY_CELL_TYPE_CLONECALL_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p35plot <- renderPlot({
    showModal(modalDialog("Plotting figure 35..", footer=NULL))
    p <- data()[['p35plots']]
    removeModal()
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p35plot, {
    screenshot(id="p35plot", filename = paste("35SCA_scImmune_CLONOTYPES_BY_CELL_TYPE_GROUP_DATA_TYPE_CLONECALL_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p36plot <- renderPlot({
    showModal(modalDialog("Plotting figure 36..", footer=NULL))
    p <- data()[['p36plots']][[input$p36id]]
    removeModal()
    return(p)
  }, height = 450, width = 900)
  
  observeEvent(input$p36plot, {
    screenshot(id="p36plot", filename = paste("36SCA_scImmune_CLONO_DIVERISITY_BY_CELL_TYPE_CLONECALL_",input$p36id, sep = ""), scale = 2)
  })
  
  output$p37plot <- renderPlot({
    showModal(modalDialog("Plotting figure 37..", footer=NULL))
    p <- data()[['p37plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p37plot, {
    screenshot(id="p37plot", filename = paste("37SCA_scImmune_CLONAL_OVERLAY_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p38plot <- renderPlot({
    showModal(modalDialog("Plotting figure 38..", footer=NULL))
    p <- data()[['p38plots']]
    removeModal()
    return(p)
  }, height = 700, width = 1400)
  
  observeEvent(input$p38plot, {
    screenshot(id="p38plot", filename = paste("38SCA_scImmune_CLONAL_OVERLAY_CELL_TYPES_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p39plot <- renderPlot({
    showModal(modalDialog("Plotting figure 39..", footer=NULL))
    p <- data()[['p39plots']]
    removeModal()
    return(p)
  }, height = 900, width = 1400)
  
  observeEvent(input$p39plot, {
    screenshot(id="p39plot", filename = paste("39SCA_scImmune_CLONAL_NETWORK_CLUSTERS_",data()[['project_name']], sep = ""), scale = 2)
  })
  
  output$p40plot <- renderPlot({
    showModal(modalDialog("Plotting figure 40..", footer=NULL))
    p <- data()[['p40plots']]
    removeModal()
    return(p)
  }, height = 1000, width = 1400)
  
  observeEvent(input$p40plot, {
    screenshot(id="p40plot", filename = paste("40SCA_scImmune_CLONAL_NETWORK_CELL_TYPES_",data()[['project_name']], sep = ""), scale = 2)
  })
  
})


