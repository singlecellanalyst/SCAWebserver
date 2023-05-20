packages <- c("shiny","BiocManager","devtools","remotes","ggplot2","SDMTools","shiny","shinyjs","robustbase","rhandsontable","gridExtra","reshape2","ggridges","patchwork","gplots","plot3D","umap","viridis","cowplot","gridGraphics","ggthemes","ggrepel","scales","rJava","cytolib","flowCore","ggcyto","flowStats","openCyto","flowViz","limma","Biobase","cytofWorkflow","FlowSOM","ConsensusClusterPlus","flowWorkspace","flowDensity","vite","FLOWMAP","CytoExploreRData","CytoExploreR"); lapply(packages, library, character.only = TRUE)

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