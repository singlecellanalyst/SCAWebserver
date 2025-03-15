library(RISmed)

ckeywords <- NULL
ckeywords$scrnaseq <- c("\"scRNAseq\" OR \"scRNA-seq\" OR \"sc-RNA-seq\" OR \"(sc)RNA-seq\" OR \"single-cell RNA sequencing\" OR \"single-cell transcriptomics\" OR \"single-cell transcriptome\" OR \"single-cell gene expression\" OR \"single-cell RNA-seq\"")
ckeywords$scatac <- c("\"scATACseq\" OR \"scATAC-seq\" OR \"sc-ATAC-seq\" OR \"(sc)ATAC-seq\" OR \"single-cell ATAC sequencing\" OR \"single-cell chromatin accessibility\" OR \"single-cell ATAC-seq\" OR \"single-cell ATAC-seq\"")
ckeywords$scimmune <- c("\"scImmune profiling\" OR \"scImmune-profiling\" OR \"single cell VDJ profiling\" OR \"single cell VDJ\" OR \"single-cell VDJ\" OR \"single-cell immune profiling\"")
ckeywords$sccnv <- c("\"scCNV\" OR \"single-cell copy number variation\" OR \"single cell copy number aberrations\" OR \"single cell copy number\"")
ckeywords$spatial <- c("\"visium spatial transcriptomics\" OR \"Spatial Transcriptomics\" OR \"Spatial Transcriptome\" OR \"spatial transcriptomic\" OR \"spatial rna seq\"")
ckeywords$flow <- c("flow cytometry")
ckeywords$cytof <- c("\"mass cytometry\" OR \"cytof\" OR \"mass spectrometry\"")
ckeywords$multiomics <- c("\"single-cell Multi-omics\" OR \"single-cell multimodal\" OR \"single cell Multi-omics\" OR \"single cell Multimodal\"")

query_journals <- NULL
query_years <- NULL
query_abstracts <- NULL

for(i in 1:length(ckeywords)){
  for(j in 2008:2022){
    print(paste(i,j))
    ngs_search <- NULL
    abstracts <- NULL
    fetch <- NULL
    ngs_search <- EUtilsSummary(ckeywords[[i]],type = "esearch", db = "pubmed",datetype = "pdat",
                                retmax = 30000,mindate=j, maxdate=j)
    if(ngs_search@count != 0){
      fetch <- EUtilsGet(ngs_search, type = "efetch", db = "pubmed")
      abstracts <- data.frame(omics = names(ckeywords)[i],
                              year = fetch@YearPubmed,
                              title = fetch@ArticleTitle,
                              abstract = fetch@AbstractText, 
                              journal = fetch@ISOAbbreviation,
                              DOI = fetch@PMID)
      query_abstracts <- rbind(query_abstracts, abstracts)
      ngs_records <- NULL
      ngs_records <- EUtilsGet(ngs_search)
      years <- NULL
      years <- YearPubmed(ngs_records)
      ngs_pubs_count <- NULL
      ngs_pubs_count <- as.data.frame(table(years))
      query_years <- rbind(query_years, data.frame(omics = names(ckeywords)[i], ngs_pubs_count))
      journal <- NULL
      journal <- ngs_records@ISOAbbreviation
      query_journals <- rbind(query_journals, data.frame(Year = i, Omics = names(ckeywords)[i], journal))
      Sys.sleep(5)
    }
  }
  
  saveRDS(query_journals, paste("query_journals_tools_",names(ckeywords)[i],".RDS", sep = ""))
  saveRDS(query_years, paste("query_years_tools_",names(ckeywords)[i],".RDS", sep = ""))
  saveRDS(query_abstracts, paste("query_abstracts_tools_",names(ckeywords)[i],".RDS", sep = ""))
  
}

