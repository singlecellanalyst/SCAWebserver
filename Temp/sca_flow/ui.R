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
    
    title = "Flow | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>Flow Cytometry Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "Flow Cytometry Online Analysis Platform",
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
                        'Upload a zip file containing Flow Cytometry FCS files (please click submit only after all your files have been uploaded):'),
                    multiple = TRUE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                downloadButton("downloadExample1", label = "Flow Input Example"),
                # downloadButton("example2", label = "gene-cell.csv example"),
                em("Alternatively, right click on the button and choose 'Save Link As'"),
                a(
                  "(dataset taken from: Dillon Hammill,2021).",
                  href = "https://github.com/DillonHammill/CytoExploreRData/tree/master/inst/extdata/Activation",
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
                downloadButton("downloadExample2", label = "Flow Metadata Example"),
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
                sidebarLayout(
                    sidebarPanel(
                        h3("ARCSINH DENSITY OF CHANNELS"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p1id",
                            label = "Choose a channel to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p1plot", "Download")
                    ),
                    mainPanel(plotOutput("p1plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel( # width = 12,
                                 h3("GATING STRATEGY PART 1"),
                                 selectInput(
                                     "p21id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 ),
                                 actionButton("p21plot", "Download")
                    ),
                    mainPanel(plotOutput("p21plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                                 h3("GATING STRATEGY PART 2"),
                                 selectInput(
                                     "p22id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 ),
                                 actionButton("p22plot", "Download")
                    ),
                    mainPanel(plotOutput("p22plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                                 h3("GATING STRATEGY PART 3"),
                                 selectInput(
                                     "p23id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 ),
                                 actionButton("p23plot", "Download")
                    ),
                    mainPanel(plotOutput("p23plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                                 h3("GATING STRATEGY PART 4"),
                                 selectInput(
                                     "p24id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 ),
                                 actionButton("p24plot", "Download")
                    ),
                    mainPanel(plotOutput("p24plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                h3("GATING SUMMARY TABLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                fluidRow(column(12, DTOutput('p2table'))),
                br(),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("POST-GATING ARCSINH DENSITY OF CHANNELS"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p3id",
                            label = "Choose a channel to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p3plot", "Download")
                    ),
                    mainPanel(plotOutput("p3plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                # h3("CELL COUNT AT EVERY GATING STAGE",
                #    style={'background-color:#EDF0F1;
                #        text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                # plotOutput("p4plot", width = "1200px", height = "800px"),
                # actionButton("p4plot", "Screenshot"),
                sidebarLayout(
                  sidebarPanel(
                    h3("CELL COUNT AT EVERY GATING STAGE"),
                    # helpText("Change Variables to Display"),
                    selectInput(
                      "p4id",
                      label = "Choose a category to display",
                      choices = "",
                      selected = ""
                    ),
                    actionButton("p4plot", "Download")
                  ),
                  mainPanel(plotOutput("p4plot", width = "900px", height = "500px")),
                  #plotlyOutput
                  position = "right",
                ),
                
                br(),
                br()
            ),
            tabPanel(
                "2.Diagnostic Plots",
                value = "section2",
                h2("Diagnostic Plots"),
                h5("Diagnostic plots showing sample distance and marker information. Channels which have no markers were removed, the remaining proceed on with the downstream analysis."),
                br(),
                h3("SAMPLE DISTANCE (MDS)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p5plot", width = "1200px", height = "800px"),
                br(),
                actionButton("p5plot", "Screenshot"),
                br(),
                br(),
                h3("SAMPLE DISTANCE (PCA)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p6plot", width = "1200px", height = "800px"),
                actionButton("p6plot", "Screenshot"),
                br(),
                br(),
                h3("ARCSINH MEDIAN EXPRESSION OF MARKERS FOR EACH SAMPLE\n(Z-Score scaling for both rows and columns)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p7plot", width = "1200px", height = "600px"),
                actionButton("p7plot", "Screenshot"),
                br(),
                br(),
                h3("DENDROGRAM BASED ON ARCSINH MEDIAN EXPRESSION\n(Z-Score scaling for both rows and columns)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p8plot", width = "1000px", height = "1000px"),
                actionButton("p8plot", "Screenshot"),
                br(),
                br(),
                
                sidebarLayout(
                  sidebarPanel(
                    h3("MARKER RANKING"),
                    # helpText("Change Variables to Display"),
                    selectInput(
                      "p9id",
                      label = "Choose a category to display",
                      choices = "",
                      selected = ""
                    ),
                    actionButton("p9plot", "Download")
                  ),
                  mainPanel(plotOutput("p9plot", width = "900px", height = "600px")),
                  #plotlyOutput
                  position = "right",
                ),
                br(),
                br()
            ),
            tabPanel(
                "3.Dimension Reduction & Clustering",
                value = "section3",
                h2("Dimension Reduction and Clustering"),
                h5("Visualisation of dimension reduction and clustering results."),
                br(),
                h3("UMAP CLUSTERS OF ALL SAMPLES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p10plot", width = "1200px", height = "800px"),
                br(),
                actionButton("p10plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP CLUSTERS - SEPARATE BY SAMPLES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p11plot", width = "1200px", height = "600px"),
                actionButton("p11plot", "Screenshot"),
                br(),
                br(),
                h3("PCA CLUSTERS OF ALL SAMPLES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p12plot", width = "1200px", height = "800px"),
                actionButton("p12plot", "Screenshot"),
                br(),
                br(),
                h3("ARCSINH MEDIAN EXPRESSION OF MARKERS IN EACH CLUSTER\n(Z-Score scaling for both rows and columns)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p13plot", width = "800px", height = "1000px"),
                actionButton("p13plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - ARCSINH EXPRESSION OF MARKERS IN ALL CELLS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p14plot", width = "1400px", height = "1000px"),
                actionButton("p14plot", "Screenshot"),
                br(),
                br(),
                h3("SAMPLE PROPORTION IN EACH CLUSTER",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p15plot", width = "1200px", height = "1200px"),
                actionButton("p15plot", "Screenshot"),
                br(),
                br()
            ),
            tabPanel(
                "4.Marker Density Plots (Clusters)",
                value = "section4",
                h2("Marker Density Plots (Clusters)"),
                h5("Marker expression density plots in each cluster. Density plot in red represents the density distribution of a particular marker in a particular cluster whereas density plot in black represents its density distribution in the entire dataset (all clusters). Markers were ranked in decreasing order based on the fold-change of the median expression of the marker in the cluster compared to the median expression of the marker in the entire dataset. Vertical lines represent the median expression value of the respective distributions. Here, only the markers whose expression is positive (arcsinh expression >0) in the selected cluster and the number of cells in the cluster having to positively express these markers are more than 3 cells, are plotted."),
                sidebarLayout(
                    sidebarPanel(
                        h3("ARCSINH DENSITY OF MARKERS IN EACH CLUSTER"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p16id",
                            label = "Choose a cluster to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p16plot", "Download")
                    ),
                    mainPanel(plotOutput("p16plot", width = "900px", height = "500px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
            )
        )
    )
))
