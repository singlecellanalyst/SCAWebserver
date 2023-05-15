library(shiny)
library(shinythemes)
library(ggplot2)
library(DT)
library(shiny)
library(shinyWidgets)
library(shinyscreenshot)

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
    
    title = "CyTOF | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>CyTOF Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "CyTOF Online Analysis Platform",
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
                        'Upload a zip file containing CyTOF CD45+ FCS files (please click submit only after all your files have been uploaded):'),
                    multiple = TRUE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                downloadButton("downloadExample1", label = "CyTOF Input Example"),
                # downloadButton("example2", label = "gene-cell.csv example"),
                em("Alternatively, right click on the button and choose 'Save Link As'"),
                a(
                  "(dataset taken from: Nowicka, M., et al. (2017)).",
                  href = "https://doi.org/10.12688/f1000research.11622.3",
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
                downloadButton("downloadExample2", label = "CyTOF Metadata Example"),
                em("Alternatively, right click on the button and save link as."),
                br(),
                br(),
                checkboxInput(
                    "termscheck",
                    a(
                      "I have read and agreed to the Terms and Conditions of the Single Cell Analyst platform.",
                      href = "https://www.singlecellanalyst.org/termsandconditions",
                      target="_blank",
                      style = "color:black"
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
                "1.Basic Metrics",
                value = "section1",
                h2("Basic Data Metrics"),
                h5("Plots showing basic input data metrics."),
                br(),
                h3("CELL COUNTS PER SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p1plot", width = "1400px", height = "1000px"),
                br(),
                actionButton("p1plot", "Screenshot"),
                br(),
                br(),
                h3("ARCSINH MARKER EXPRESSION FOR EACH GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p2plot", width = "1400px", height = "1000px"),
                br(),
                actionButton("p2plot", "Screenshot"),
                br(),
                br(),
                h3("MDS - SAMPLE DISTANCE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p3plot", width = "1400px", height = "800px"),
                br(),
                actionButton("p3plot", "Screenshot"),
                br(),
                br(),
                h3("PCA - SAMPLE DISTANCE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p4plot", width = "1400px", height = "800px"),
                br(),
                actionButton("p4plot", "Screenshot"),
                br(),
                br()
            ),
            tabPanel(
                "2.Marker Phenotypes",
                value = "section2",
                h2("Marker Phenotypes"),
                h5("Phenotypic discoveries of markers in each sample or each group. Non-marker channels are removed, and only markers shown in the plots are considered in the downstream analyses."),
                br(),
                h3("MARKER ARCSINH MEDIAN EXPRESSION - SAMPLES\n(Z-Score scaling for both rows and columns)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p5plot", width = "1400px", height = "800px"),
                br(),
                actionButton("p5plot", "Screenshot"),
                br(),
                br(),
                h3("BOX PLOTS OF MARKERS FOR EACH GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p6plot", width = "1400px", height = "1400px"),
                actionButton("p6plot", "Screenshot"),
                br(),
                br(),
                h3("DENDROGRAM BASED ON ARCSINH MEDIAN EXPRESSION\n(Z-Score scaling for both rows and columns)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p7plot", width = "1000px", height = "1000px"),
                actionButton("p7plot", "Screenshot"),
                br(),
                br(),
                h3("NON-REDUNDANCY SCORES - MARKER IMPORTANCE RANKING",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p8plot", width = "1400px", height = "700px"),
                actionButton("p8plot", "Screenshot"),
                br(),
                br(),
                h3("CHOSEN K NUMBER OF CLUSTERS BASED ON PAC VALUE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p9plot", width = "1400px", height = "700px"),
                actionButton("p9plot", "Screenshot"),
                br(),
                br()
            ),
            tabPanel(
                "3.Dimension Reduction & Clustering",
                value = "section3",
                h2("Dimension Reduction and Clustering"),
                h5("Visualisation of dimension reduction and clustering results.Non-marker channels are removed, and only markers shown in the plots are considered in the downstream analyses."),
                br(),
                h3("ARCSINH MEDIAN EXPRESSION OF MARKERS IN EACH CLUSTER\n(Z-Score scaling for both rows and columns)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p10plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p10plot", "Screenshot"),
                br(),
                br(),
                h3("BOXPLOT - ARCSINH MEDIAN EXPRESSION OF MARKERS PER CLUSTER",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p11plot", width = "1400px", height = "1400px"),
                actionButton("p11plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p12plot", width = "1400px", height = "800px"),
                actionButton("p12plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    h3("DENSITY OF MARKERS IN EACH CLUSTER AND GROUP"),
                    # helpText("Change Variables to Display"),
                    selectInput(
                      "p11bid",
                      label = "Choose a cluster to display",
                      choices = "",
                      selected = ""
                    ),
                    actionButton("p11bplot", "Download")
                  ),
                  mainPanel(plotOutput("p11bplot", width = "1200px", height = "900px")),
                  #plotlyOutput
                  position = "right",
                ),
                br(),
                br(),
                h3("UMAP - BY SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p13plot", width = "1400px", height = "800px"),
                actionButton("p13plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - BY CLUSTER",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p14plot", width = "1400px", height = "800px"),
                actionButton("p14plot", "Screenshot"),
                br(),
                br(),
                h3("PCA - BY CLUSTER",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p15plot", width = "1400px", height = "800px"),
                actionButton("p15plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("ARCSINH DENSITY OF MARKERS IN EACH GROUP"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p16id",
                            label = "Choose a group to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p16plot", "Download")
                    ),
                    mainPanel(plotOutput("p16plot", width = "1000px", height = "1400px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                h3("PROPORTION OF SAMPLES IN EACH CLUSTER (PIE)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p17plot", width = "1400px", height = "1400px"),
                actionButton("p17plot", "Screenshot"),
                br(),
                br(),
                h3("PROPORTION OF SAMPLES IN EACH CLUSTER (STACKED BAR PLOT)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p17bplot", width = "1400px", height = "800px"),
                actionButton("p17bplot", "Screenshot"),
                br(),
                br(),
                h3("PROPORTION OF CLUSTERS IN EACH SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p18plot", width = "1400px", height = "800px"),
                actionButton("p18plot", "Screenshot"),
                br(),
                br(),
                h3("SOM NODES - tSNE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p19plot", width = "1400px", height = "800px"),
                actionButton("p19plot", "Screenshot"),
                br(),
                br(),
                h3("SOM NODES PHENOTYPES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                fluidRow(column(12, DTOutput('p20table'))),
                br(),
                br(),
                h3("SOM NODES - UMAP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p21plot", width = "1400px", height = "800px"),
                actionButton("p21plot", "Screenshot"),
                br(),
                br(),
                h3("SOM NODES - PCA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p22plot", width = "1400px", height = "800px"),
                actionButton("p22plot", "Screenshot"),
                br(),
                br(),
            ),
            tabPanel(
                "4.Marker Density Plots (Clusters)",
                value = "section4",
                h2("Marker Density Plots (Clusters)"),
                h5("Marker expression density plots in each cluster. Density plot in red represents the density distribution of a particular marker in a particular cluster whereas density plot in black represents its density distribution in the entire dataset (all clusters). Markers were ranked in decreasing order based on the fold-change of the median expression of the marker in the cluster compared to the median expression of the marker in the entire dataset. Vertical lines represent the median expression value of the respective distributions. Here, only the markers whose expression is positive (arcsinh expression >0) in the selected cluster and the number of cells in the cluster having to positively express these markers are more than 3 cells, are plotted."),


                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("DENSITY OF MARKERS IN EACH CLUSTER"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p23id",
                            label = "Choose a cluster to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p23plot", "Download")
                    ),
                    mainPanel(plotOutput("p23plot", width = "900px", height = "1200px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br()
            )
        )
    )
))
