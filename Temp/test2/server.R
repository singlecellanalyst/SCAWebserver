packages <- c("remotes","devtools","BiocManager","R.utils","patchwork","AnnotationFilter","BiocGenerics","GenomicFeatures","GenomicRanges","IRanges","Rsamtools","S4Vectors","ggbio","motifmatchr","AnnotationDbi","Seurat","GenomeInfoDb","JASPAR2020","TFBSTools","EnsDb.Hsapiens.v86","BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19","Signac","shiny","patchwork","hdf5r","ComplexHeatmap","rJava"); lapply(packages, library, character.only = TRUE)

function(input, output) {
  
  output$main_plot <- renderPlot({
    
    hist(faithful$eruptions,
         probability = TRUE,
         breaks = as.numeric(input$n_breaks),
         xlab = "Duration (minutes)",
         main = "Geyser eruption duration")
    
    if (input$individual_obs) {
      rug(faithful$eruptions)
    }
    
    if (input$density) {
      dens <- density(faithful$eruptions,
                      adjust = input$bw_adjust)
      lines(dens, col = "blue")
    }
    
  })
}