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
          }, 1500000000)
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

    title = "scATAC | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>Single Cell ATAC Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "Single Cell ATAC Online Analysis Platform",
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
                        'Upload input a zipped file containing 10X scATAC 1) peak_bc_matrix.h5, 2) fragments.tsv.gz and 3) its index .tbi file, 4) per_barcode_metrics.csv or singlecell.csv files. (Please submit after the upload is completed)'),
                    multiple = FALSE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                # downloadButton("downloadExample1", label = "scATAC Input Example"),
                # downloadButton("example2", label = "gene-cell.csv example"),
                actionButton("example1", label = "scATAC Input Example", onclick=paste0("window.open('https://www.dropbox.com/s/vf3e2zeckbqhsfe/SCA_scATAC_Example_From_10X.zip?dl=1','_self')")),
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
                        'Upload a meta file in .csv/.txt format for your data (format see example, please submit after the upload is completed)'
                    ),
                    accept = c(
                        'text/csv',
                        'text/txt',
                        'text/comma-separated-values,text/plain',
                        '.csv',
                        '.txt'
                    )
                ),
                downloadButton("downloadExample2", label = "scATAC Metadata Example"),
                em("Alternatively, right click on the button and save link as."),
                br(),
                br(),
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
                "1.Basic Metrics",
                value = "section1",
                h2("Quality Control Metrics"),
                h5("Plots showing different QC metrics and how input data looks like before and after filtering has been done."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TSS ENRICHMENT SCORES"),
                        selectInput(
                            "p1id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p1plot", "Download")
                    ),
                    mainPanel(plotOutput("p1plot", height = "450px", width = "900px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("FRAGMENT LENGTH PERIODICITY"),
                        selectInput(
                            "p2id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p2plot", "Download")
                    ),
                    mainPanel(plotOutput("p2plot", height = "450px", width = "900px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("PEAKS QC"),
                        selectInput(
                            "p3id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p3plot", "Download")
                    ),
                    mainPanel(plotOutput("p3plot", height = "450px", width = "900px")),
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "2.Dimension Reduction & Clustering",
                value = "section2",
                h2("Dimension Reduction & Clustering"),
                h5("Plots showing dimension reduction and clustering steps for each scATAC sample."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("CORRELATION BETWEEN DEPTH VS LSI"),
                        selectInput(
                            "p4id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p4plot", "Download")
                    ),
                    mainPanel(plotOutput("p4plot", width = "900px", height = "450px")),

                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP - CLUSTERS"),
                        selectInput(
                            "p5id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p5plot", "Download")
                    ),
                    mainPanel(plotOutput("p5plot", width = "900px", height = "600px")),

                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("TOP PEAKS IN EACH CLUSTER"),
                                 selectInput(
                                     "p6id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p6table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("CLOSEST GENE AT TOP PEAKS"),
                                 selectInput(
                                     "p8id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p8table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("TOP GENE ACTIVITY"),
                                 selectInput(
                                     "p9id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p9table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "3.Phenotypic Discovery",
                value = "section3",
                h2("Phenotypic Discovery"),
                h5("Insertion coverage, footprinting and motif analyses for top peaks in clusters."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("Tn5 INSERTION COVERAGE - TOP PEAK + ACTIVITY"),
                        selectInput(
                            "p10id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput(
                            "p11id",
                            label = "Choose a cluster to display",
                            value = 0,min = 0,max = 10,step = 1
                        ),
                        actionButton("p10plot", "Download")
                    ),
                    mainPanel(plotOutput("p10plot", width = "900px", height = "700px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("FOOT PRINTING ANALYSIS"),
                        selectInput(
                            "p12id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p12plot", "Download")
                    ),
                    mainPanel(plotOutput("p12plot", width = "900px", height = "900px")),

                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP ENRICHED MOTIFS"),
                        selectInput(
                            "p14id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p14plot", "Download")
                    ),
                    mainPanel(plotOutput("p14plot", width = "900px", height = "700px")),

                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("ENRICHED MOTIFS"),
                                 selectInput(
                                     "p13id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p13table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
                # sidebarLayout(
                #     sidebarPanel(
                #         h3("UMAP WITH TOP ENRICHED MOTIFS"),
                #         selectInput(
                #             "p15id",
                #             label = "Choose a sample to display",
                #             choices = "",
                #             selected = ""
                #         ),
                #         actionButton("p15plot", "Download")
                #     ),
                #     mainPanel(plotOutput("p15plot", width = "900px", height = "1400px")),
                #     position = "right",
                # ),
                # br(),
                # br(),
            ),
            tabPanel(
                "4.Merged Analysis",
                value = "section4",
                h2("Merged Analysis"),
                h5("Merged analysis for project with more than one sample."),
                br(),
                h3("UMAP - SAMPLES & CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: left; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p16plot", width = "1400px", height = "700px"),
                actionButton("p16plot", "Screenshot"),
                br(),
                br(),
                h3("TOP COMMON PEAKS IN CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: left; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p17plot", width = "1200px", height = "1600px"),
                actionButton("p17plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP PEAK IN EACH CLUSTER"),
                        # selectInput(
                        #     "p18id",
                        #     label = "Choose a cluster to display",
                        #     choices = "",
                        #     selected = ""
                        # ),
                        sliderInput(
                            "p18id",
                            label = "Choose a cluster to display",
                            value = 0,min = 0,max = 10,step = 1
                        ),
                        actionButton("p18plot", "Download")
                    ),
                    mainPanel(plotOutput("p18plot", width = "900px", height = "700px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("Tn5 INSERTION COVERAGE - TOP PEAK + ACTIVITY"),
                        # selectInput(
                        #     "p19id",
                        #     label = "Choose a cluster to display",
                        #     choices = "",
                        #     selected = ""
                        # ),
                        sliderInput(
                            "p19id",
                            label = "Choose a cluster to display",
                            value = 0,min = 0,max = 10,step = 1
                        ),
                        actionButton("p19plot", "Download")
                    ),
                    mainPanel(plotOutput("p19plot", width = "900px", height = "800px")),

                    position = "right",
                ),
                br(),
                br()
            )
        )
    )
))
