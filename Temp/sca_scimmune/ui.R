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
    
    title = "scImmune | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>scImmune Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "scImmune Online Analysis Platform",
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
                        'Upload input a zipped file containing 10X Immune Profiling 1) filtered contig, 2) filtered consensus, 3) filtered gene matrix feature_bc_matrix.tar.gz or .h5 files (please click submit only after all your files have been uploaded):'),
                    multiple = TRUE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                # downloadButton("downloadExample1", label = "scImmune Input Example"),
                # downloadButton("example2", label = "gene-cell.csv example"),
                actionButton("example1", label = "scImmune Input Example", onclick=paste0("window.open('https://www.dropbox.com/s/il954lexu6hn6fe/SCA_scImmune_Example_From_10X.zip?dl=1','_self')")),
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
                # downloadButton("downloadExample2", label = "scImmune Metadata Example"),
                actionButton("example2", label = "scImmune Metadata Example", onclick=paste0("window.open('https://www.dropbox.com/s/0lbxq9puegb714i/SCA_scImmune_Metadata_Example.csv?dl=1','_self')")),
                em("Alternatively, right click on the button and save link as."),
                br(),
                br(),
                selectInput(
                    "integration_method",
                    "4. Choose an integration method (only applies when multiple samples are uploaded):",
                    choices = c("Seurat","Harmony"),
                    selected = "Seurat"),
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
                h2("Basic Data Metrics"),
                h5("Plots showing basic input data metrics."),
                br(),
                h3("NUMBER OF UNIQUE CLONOTYPES (SAMPLES & GROUPS)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p1plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p1plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("CLONOTYPE ABUNDANCE"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p18id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p18plot", "Download")
                    ),
                    mainPanel(plotOutput("p18plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                h3("DISTRIBUTION OF CDR3 LENGTHS (AA & NT)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p2plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p2plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("CONTIG LENGTH"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p20id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p20plot", "Download")
                    ),
                    mainPanel(plotOutput("p20plot", width = "900px", height = "600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("AA ABUNDANCE"),
                                 selectInput(
                                     "p19id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p19table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("REPERTOIRE SPECTRATYPES"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p13id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p13plot", "Download")
                    ),
                    mainPanel(plotOutput("p13plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "2.Clonal Analysis",
                value = "section2",
                h2("Clonal Analysis"),
                h5("Plots showing clonal proportions and clonal space homeostatis."),
                br(),
                h3("CLONAL SPACE HOMEOSTASIS - CDR3 AA REGION (SAMPLES & GROUPS)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p6plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p6plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("CLONAL SPACE HOMEOSTASIS & PROPORTIONS"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p23id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p23plot", "Download")
                    ),
                    mainPanel(plotOutput("p23plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                h3("RARE CLONAL PROPORTIONS (SAMPLES & GROUPS)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p5plot", width = "1400px", height = "700px"),
                actionButton("p5plot", "Screenshot"),
                br(),
                br(),
            ),
            tabPanel(
                "3.Repertoire Discovery",
                value = "section3",
                h2("Repertoire Discovery"),
                h5("Visualisation of TCR/BCR repertoires."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP CLONOTYPES (SAMPLEWISE)"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p15id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p15plot", "Download")
                    ),
                    mainPanel(plotOutput("p15plot", width = "900px", height = "1600px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("OVERALL TOP CLONOTYPES"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p21id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p21plot", "Download")
                    ),
                    mainPanel(plotOutput("p21plot", width = "900px", height = "800px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                h3("REPERTOIRE OVERLAPS (VDJ GENES + CDR3 NUCLEOTIDES)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p7plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p7plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("REPERTOIRE OVERLAPS"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p24id",
                            label = "Choose a cell type to display",
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
                h3("HEATMAP - REPERTOIRE OVERLAPS (CDR3 AA)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p17plot", width = "1200px", height = "900px"),
                actionButton("p17plot", "Screenshot"),
                br(),
                br(),
                h3("CIRCOS - REPERTOIRE OVERLAPS (CDR3 AA)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p16plot", width = "1200px", height = "1200px"),
                actionButton("p16plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("HIERARCHICAL CLUSTERING (SAMPLES)"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p25id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p25plot", "Download")
                    ),
                    mainPanel(plotOutput("p25plot", width = "900px", height = "500px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "4.Gene Usage",
                value = "section4",
                h2("Gene Usage"),
                h5("Visualisation of VDJ(C) gene usage."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP GENES USAGE"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p9id",
                            label = "Choose a gene to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p9plot", "Download")
                    ),
                    mainPanel(plotOutput("p9plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("V GENE USAGE"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p22id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p22plot", "Download")
                    ),
                    mainPanel(plotOutput("p22plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("GENE USAGE CORRELATION"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p10id",
                            label = "Choose a gene to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p10plot", "Download")
                    ),
                    mainPanel(plotOutput("p10plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("SAMPLE DISTANCE"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p12id",
                            label = "Choose a gene to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p12plot", "Download")
                    ),
                    mainPanel(plotOutput("p12plot", width = "900px", height = "400px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("HIERARCHICAL CLUSTERING"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p11id",
                            label = "Choose a gene to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p11plot", "Download")
                    ),
                    mainPanel(plotOutput("p11plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "5.Diversity Estimation",
                value = "section5",
                h2("Diversity Estimation"),
                h5("Estimate clonal diversity using different methods."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("CLONAL DIVERSITY"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p26id",
                            label = "Choose a cell type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p26plot", "Download")
                    ),
                    mainPanel(plotOutput("p26plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("DIVERSITY ESTIMATION ANALYSIS"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p14id",
                            label = "Choose a method",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p14plot", "Download")
                    ),
                    mainPanel(plotOutput("p14plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP K-MERS"),
                        # helpText("Change Variables to Display"),
                        sliderInput(
                            "kmerid",label = "Choose a k to show its top k-mers",
                            min = 2, max = 15, value = 2, step = 1),
                        actionButton("kmerplot", "Download")
                    ),
                    mainPanel(plotOutput("kmerplot", width = "900px", height = "700px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("K-MER SEQUENCE MOTIFS"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "sampleid",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput(
                            "kmer2id",label = "Choose a k to show its top k-mers",
                            min = 2, max = 15, value = 2, step = 1),
                        actionButton("kmer2plot", "Download")
                    ),
                    mainPanel(plotOutput("kmer2plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "6.Integrative Phenotypes",
                value = "section6",
                h2("Integrative Phenotypes"),
                h5("Integration of scImmune data with scRNA-Sequencing data."),
                br(),
                h3("UMAP - CELL TYPE PREDICTION",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p31plot", width = "1400px", height = "1000px"),
                actionButton("p31plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - CELL TYPE BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p32plot", width = "1400px", height = "700px"),
                actionButton("p32plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - AUTO CLUSTERING",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p29plot", width = "1400px", height = "1000px"),
                actionButton("p29plot", "Screenshot"),
                br(),
                br(),
                DTOutput('demarkers'),
                br(),
                br(),
                h3("UMAP - CLUSTERS BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p30plot", width = "1400px", height = "700px"),
                actionButton("p30plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - INTEGRATION BY SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p27plot", width = "1400px", height = "1000px"),
                actionButton("p27plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP/tSNE - BY TCR/BCR GROUPS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p28plot", width = "1400px", height = "700px"),
                actionButton("p28plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - CLONAL HOMEOSTASIS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p33plot", width = "1400px", height = "1000px"),
                actionButton("p33plot", "Screenshot"),
                br(),
                br(),
                h3("CLONAL HOMEOSTASIS - CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p34plot", width = "1400px", height = "700px"),
                actionButton("p34plot", "Screenshot"),
                br(),
                br(),
                h3("CLONAL HOMEOSTASIS - CLUSTERS BY GROUPS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p35plot", width = "1400px", height = "1000px"),
                actionButton("p35plot", "Screenshot"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("CLONAL DIVERSITY - CELL TYPES"),
                        # helpText("Change Variables to Display"),
                        selectInput(
                            "p36id",
                            label = "Choose a type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p36plot", "Download")
                    ),
                    mainPanel(plotOutput("p36plot", width = "900px", height = "450px")),
                    #plotlyOutput
                    position = "right",
                ),
                br(),
                br(),
                h3("CLONAL FREQUENCIES IN CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p37plot", width = "1400px", height = "700px"),
                actionButton("p37plot", "Screenshot"),
                br(),
                br(),
                h3("CLONAL FREQUENCIES IN CELL TYPES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p38plot", width = "1400px", height = "700px"),
                actionButton("p38plot", "Screenshot"),
                br(),
                br(),
                h3("CLONAL NETWORK IN CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p39plot", width = "1400px", height = "1000px"),
                actionButton("p39plot", "Screenshot"),
                br(),
                br(),
                h3("CLONAL NETWORK IN CELL TYPES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p40plot", width = "1400px", height = "1000px"),
                actionButton("p40plot", "Screenshot"),
                br(),
                br(),
            )
        )
    )
))
