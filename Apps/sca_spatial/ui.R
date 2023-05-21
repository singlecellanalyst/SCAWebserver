library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(shinyWidgets)
library(shinyscreenshot)

cc_choices <- c("No Cell Cycle Regression",
                "G2M and S Phase Scores Regression",
                "G2M and S Phase Difference Regression")

shinyUI(fluidPage(
    
    tags$head(
        HTML(
            "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 150000000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
        ),
        tags$style('body {font-size: 15px;}')
    ),
    textOutput("keepAlive"),
    
    title = "Spatial | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>Spatial Transcriptomics Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "Spatial Transcriptomics Online Analysis Platform",
        align = "center"
    ),
    h5(
        "SCA aims to facilitate researchers around the globe to utilize its one-stop single cell analysis platform and deliver breakthrough results.",
        align = "center"
    ),
    navbarPage(
        "SCA",
        id = "nav",
        # title=div(img(src="favicon/favicon-32x32.png", style = "margin:-6px 0px"), ""),
        tabPanel(
            "Upload",
            # titlePanel("Project Name"),
            wellPanel(style = "background: white",
                # sidebarPanel(
                textInput(inputId = "project_name", "1. Project Name:",
                          value = "", placeholder = "Create a project name"),
                # br(),
                strong("2. Choose Files:"),
                fileInput(
                    'files',
                    h5(
                        'Upload input a .zip file containing 10X Spatial h5 file(s) and spatial image folder(for multiple samples, place spatial files of each sample in their separate folders and zip them together):'),
                    multiple = TRUE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                downloadButton("downloadExample1", label = "Spatial Input Example"),
                # downloadButton("example2", label = "gene-cell.csv example"),
                em("Alternatively, right click on the button and choose 'Save Link As'"),
                a(
                  "(dataset taken from: 10X Genomics).",
                  href = "https://www.10xgenomics.com/resources/datasets",
                  target="_blank",
                  style = "color:#333333; font-style:italic",
                ),
                br(),
                br(),
                strong("3. Choose Meta Data:"),
                fileInput(
                    'phenodata',
                    h5(
                        'Upload a meta file in .csv/.txt format for your data (format see example)'
                    ),
                    accept = c(
                        'text/csv',
                        'text/txt',
                        'text/comma-separated-values,text/plain',
                        '.csv',
                        '.txt'
                    )
                ),
                downloadButton("downloadExample2", label = "Spatial Metadata Example"),
                em("Alternatively, right click on the button and save link as."),
                br(),
                br(),
                strong("4.Choose scRNA-Seq file"),
                fileInput(
                    'scrnaseq',
                    h5(
                        'Upload a scRNA-Sequencing 10X .h5 format or Seurat .RDS file for integration with spatial data. Current pipeline only supports single sample integration (optional)'
                    ),
                    accept = c(
                        'application/.h5',
                        '.h5','.H5',
                        'application/.RDS',
                        '.rds','.RDS'
                    )
                ),
                checkboxInput(
                    "termscheck",
                    a(
                      "I have read and agreed to the Terms and Conditions of the Single Cell Analyst platform.",
                      href = "https://www.singlecellanalyst.org/termsandconditions",
                      target="_blank",
                      style = "color:black",
                    ),
                    value = FALSE
                ),
                br(),
                fluidRow(
                    column(
                        12,
                        offset = 0,
                        align = "center",
                        id = "submit",
                        actionButton('submit', 'SUBMIT')
                    )
                )
            ),
        ),
        # br(),
        navbarMenu(
            "Analysis Results",
            menuName = "results",
            tabPanel(
                "1.QC Metrics",
                value = "qcmetrics",
                h2("Quality Control Metrics"),
                h5("Plots showing different QC metrics and how input data looks like before any filtering or normalization is done."),
                sidebarLayout(
                    sidebarPanel(
                        h3("SPATIAL QC CHECK"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p1id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p1plot", "Download")
                    ),
                    
                    mainPanel(plotOutput("p1plot", width = "900px", height = "500px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP SPATIALLY VARIABLE FEATURES"),
                        h4("(2 Different Alpha Levels)"),
                        selectInput(
                            "p2id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p2plot", "Screenshot")
                    ),
                    
                    mainPanel(plotOutput("p2plot", width = "900px", height = "500px")),
                    position = "right",
                ),
                br(),
                br()
            ),
            tabPanel(
                "2.Dimension Reduction and Clustering",
                value = "drcluster",
                h2("Dimension Reduction and Clustering"),
                h5("Unsupervised clustering of spatial samples, colored by clusters or sample origins."),
                br(),
                h3("UMAP - CLUSTERS AND SAMPLES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p3plot", width = "1400px", height = "513px"),
                actionButton("p3plot", "Screenshot"),
                br(),
                br(),
                # sidebarLayout(
                #     sidebarPanel(
                #         width = "100%",
                #         h3("UMAP: CLUSTERS AND SAMPLES"),
                #         actionButton("p3plot", "Screenshot")
                #     ),
                #     
                #     mainPanel(plotOutput("p3plot", width = "1400px", height = "650px")),
                #     position = "left",
                # ),
                h3("UMAP - CLUSTERS AND SAMPLES (SEPARATED BY SAMPLES)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p4plot", width = "1400px", height = "100%"),
                actionButton("p4plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP - CLUSTERS BY SAMPLES"),
                        selectInput(
                            "p5id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p5plot", "Screenshot")
                    ),
                    
                    mainPanel(plotOutput("p5plot", width = "900px", height = "700px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("EACH CLUSTER ON SLIDE"),
                        selectInput(
                            "p6id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p6plot", "Screenshot")
                    ),
                    
                    mainPanel(plotOutput("p6plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("SAMPLE CLUSTERS ON UMAP"),
                        selectInput(
                            "p7id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p7plot", "Screenshot")
                    ),
                    
                    mainPanel(plotOutput("p7plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "3.Cluster Phenotype Discoveries",
                value = "clusterpheno",
                h2("Cluster Phenotype Discoveries"),
                h5("Plots showing phenotypic discoveries for all clusters in sample(s)."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("SPATIALLY MAPPED CLUSTERS"),
                        selectInput(
                            "p8id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p8plot", "Screenshot")
                    ),
                    
                    mainPanel(plotOutput("p8plot", width = "900px", height = "450px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKERS IN CLUSTERS"),
                        selectInput(
                            "p9id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        selectInput(
                            "p9cluster",
                            label = "Choose a cluster to view",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p9plot", "Screenshot")
                    ),
                    mainPanel(plotOutput("p9plot", width = "900px", height = "300px")),
                    position = "right",
                ),
                br(),
                br(),

                sidebarLayout(
                    sidebarPanel(width = 12,
                        h3("TOP MARKERS IN CLUSTERS"),
                        selectInput(
                            "p10id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p10table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "4.Cell Type Discovery",
                value = "ctpheno",
                h2("Cell Type Discovery"),
                h5("Plots showing top genes expressed in each clusters and auto-cell annotation."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP GENE IN EACH CLUSTER"),
                        selectInput(
                            "p11id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p11plot", "Screenshot")
                    ),
                    mainPanel(plotOutput("p11plot", width = "900px", height = "1100px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("AUTO CELL TYPE ANNOTATION (SINGLER)"),
                        selectInput(
                            "p12id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p12plot", "Screenshot")
                    ),
                    mainPanel(plotOutput("p12plot", width = "900px", height = "1000px")),
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "5.Integrative Omics",
                value = "multiomics",
                h2("Integrative Omics"),
                h5("Plots showing integrative phenotypic discoveries of spatial data with scRNA-Seq data."),
                br(),
                h3("scRNASEQ ONLY - UMAP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p13plot", width = "1400px", height = "1000px"),
                actionButton("p13plot", "Screenshot"),
                br(),
                br(),
                h3("CELL TYPE IDENTIFICATION - UMAP & tSNE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p14plot", width = "1400px", height = "1000px"),
                actionButton("p14plot", "Screenshot"),
                br(),
                br(),
                h3("CELL TYPE PREDICTION SCORES - CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p15plot", width = "1400px", height = "1000px"),
                actionButton("p15plot", "Screenshot"),
                br(),
                br(),
                h3("PREDICTION OF SPATIAL CELL TYPE BASED ON scRNA-SEQ DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p16plot", width = "1400px", height = "1000px"),
                actionButton("p16plot", "Screenshot"),
                br(),
                br(),
)
        )
        # Post-QC Metrics (2h)
        # Sample Phenotypes (2h)
        # Sample Pathway Analysis (2h)
        # Sample Pseudotime (2h)
        # Integrative QC (1h)
        # Integrative Phenotypes (1h)
        # Integrative Pseudotime (1h)
        
    )
))
