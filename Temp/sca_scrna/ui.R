library(shiny)
library(shinythemes)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(shinyscreenshot)
library(plotly)

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
    
    title = "scRNASEQ | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>Single Cell RNA Sequencing Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "Single Cell RNA Sequencing Online Analysis Platform",
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
                        'Upload input a zipped file containing 10X scRNASEQ filtered_matrix.h5 files or gene-to-cell raw count matrices in .csv or .txt formats.'),
                    multiple = FALSE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                downloadButton("downloadExample1", label = "scRNASEQ Input Example"),
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
                downloadButton("downloadExample2", label = "scRNASEQ Metadata Example"),
                em("Alternatively, right click on the button and save link as."),
                br(),
                br(),
                selectInput(
                    "integration_method",
                    "4. Choose an integration method (only applies when multiple samples are uploaded):",
                    choices = c("SEURAT","HARMONY"),
                    selected = "SEURAT"),
                selectInput(
                    "ccregression",
                    "5. Perform cell cycle regression:",
                    choices = c("NO_REGRESSION","TWO_PHASES","PHASE_DIFFERENCE"),
                    selected = "NO_REGRESSION"),
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
                value = "section1",
                h2("Quality Control Metrics"),
                h5("Plots showing different QC metrics and how input data looks like before any filtering or normalization is done."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("Individual Feature QC Check"),
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
                        h3("Feature vs Feature QC Check"),
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
                        h3("PCA Cell Cycle Check (samples with <2000 cells will not be processed)"),
                        selectInput(
                            "p3id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p3plot", "Download")
                    ),
                    mainPanel(plotOutput("p3plot", height = "600px", width = "900px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP Cell Cycle Check"),
                        selectInput(
                            "p4id",
                            label = "Choose a sample to display (samples with <2000 cells will not be processed)",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p4plot", "Download")
                    ),
                    mainPanel(plotOutput("p4plot", height = "600px", width = "900px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP Features Check"),
                        selectInput(
                            "p5id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p5plot", "Download")
                    ),
                    mainPanel(plotOutput("p5plot", height = "800px", width = "900px")),
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "2.Post-QC Metrics",
                value = "section2",
                h2("Post Quality Control Metrics"),
                h5("Plots showing different QC metrics and how input data looks like after filtering and normalization are done."),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("Individual Feature QC Check"),
                        selectInput(
                            "p9id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p9plot", "Download")
                    ),
                    mainPanel(plotOutput("p9plot", width = "900px", height = "450px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("PCA Cell Cycle Check (samples with <2000 cells will not be processed)"),
                        selectInput(
                            "p6id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p6plot", "Download")
                    ),
                    mainPanel(plotOutput("p6plot", width = "900px", height = "600px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP Cell Cycle Check (samples with <2000 cells will not be processed)"),
                        selectInput(
                            "p7id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p7plot", "Download")
                    ),
                    mainPanel(plotOutput("p7plot", width = "900px", height = "600px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP Features Check"),
                        selectInput(
                            "p8id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p8plot", "Download")
                    ),
                    mainPanel(plotOutput("p8plot", width = "900px", height = "800px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("HVG Check"),
                        selectInput(
                            "p11id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p11plot", "Download")
                    ),
                    mainPanel(plotOutput("p11plot", width = "900px", height = "450px")),
                    
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "3.Sample Phenotypes",
                value = "section3",
                h2("Sample Phenotype Discoveries"),
                h5("Plots showing phenotypic discoveries of each sample files."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE 3D UMAP BY CELL TYPE"),
                        selectInput(
                            "p80id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        # actionButton("p80plot", "Download")
                    ),
                    mainPanel(plotlyOutput("p80plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE 3D tSNE BY CELL TYPE"),
                        selectInput(
                            "p81id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        # actionButton("p81plot", "Download")
                    ),
                    mainPanel(plotlyOutput("p81plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP CELL ANNOTATION"),
                        selectInput(
                            "p24id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p24plot", "Download")
                    ),
                    mainPanel(plotOutput("p24plot", width = "900px", height = "600px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("tSNE CELL ANNOTATION"),
                        selectInput(
                            "p25id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p25plot", "Download")
                    ),
                    mainPanel(plotOutput("p25plot", width = "900px", height = "600px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP PC FEATURES"),
                        selectInput(
                            "p12id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p12plot", "Download")
                    ),
                    mainPanel(plotOutput("p12plot", width = "900px", height = "450px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("VISUALIZE PC1 VS PC2"),
                        selectInput(
                            "p13id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p13plot", "Download")
                    ),
                    mainPanel(plotOutput("p13plot", width = "900px", height = "600px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE 3D UMAP BY CLUSTER"),
                        selectInput(
                            "p82id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p82plot", "Download")
                    ),
                    mainPanel(plotlyOutput("p82plot", width = "900px", height = "600px")),
                    
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE 3D tSNE BY CLUSTER"),
                        selectInput(
                            "p83id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p83plot", "Download")
                    ),
                    mainPanel(plotlyOutput("p83plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("SAMPLE UMAP AUTOCLUSTERING"),
                        selectInput(
                            "p15id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p15plot", "Download")
                    ),
                    mainPanel(plotOutput("p15plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("SAMPLE tSNE AUTOCLUSTERING"),
                        selectInput(
                            "p16id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p16plot", "Download")
                    ),
                    mainPanel(plotOutput("p16plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("MEDIAN EXPRESSION OF TOP MARKERS IN CLUSTERS"),
                        selectInput(
                            "p18id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p18plot", "Download")
                    ),
                    mainPanel(plotOutput("p18plot", width = "900px", height = "800px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("TOP MARKERS IN CLUSTERS"),
                                 selectInput(
                                     "p17id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p17table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKERS IN CLUSTERS + tSNE VISUALIZATION"),
                        selectInput(
                            "p19id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p19plot", "Download")
                    ),
                    mainPanel(plotOutput("p19plot", width = "900px", height = "450px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKERS IN CLUSTERS + UMAP VISUALIZATION"),
                        selectInput(
                            "p193id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p193plot", "Download")
                    ),
                    mainPanel(plotOutput("p193plot", width = "900px", height = "450px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKER EXPRESSIONS IN EACH CLUSTER"),
                        selectInput(
                            "p20id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p20plot", "Download")
                    ),
                    mainPanel(plotOutput("p20plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKER EXPRESSIONS IN EACH CLUSTER VIOLIN PLOTS"),
                        selectInput(
                            "p22id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput(
                            "p222id",
                            label = "Choose a cluster to display",
                            min = 0, max = 10, value = 0, step = 1),
                        actionButton("p22plot", "Download")
                    ),
                    mainPanel(plotOutput("p22plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("UMAP: TOP GENE IN EACH CLUSTER"),
                        selectInput(
                            "p21id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p21plot", "Download")
                    ),
                    mainPanel(plotOutput("p21plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("VISUALIZE PC COMPONENTS"),
                        selectInput(
                            "p14id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p14plot", "Download")
                    ),
                    mainPanel(plotOutput("p14plot", width = "900px", height = "1000px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("HEATMAP: TOP GENES IN EACH CLUSTER"),
                        selectInput(
                            "p23id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p23plot", "Download")
                    ),
                    mainPanel(plotOutput("p23plot", width = "900px", height = "900px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("MEDIAN EXPRESSION TOP GENES IN ANNOTATED CELLTYPES"),
                        selectInput(
                            "p36id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p36plot", "Download")
                    ),
                    mainPanel(plotOutput("p36plot", width = "900px", height = "1200px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKERS IN CELL TYPES + tSNE VISUALIZATION"),
                        selectInput(
                            "p37id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p37plot", "Download")
                    ),
                    mainPanel(plotOutput("p37plot", width = "900px", height = "450px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKERS IN CELL TYPES + UMAP VISUALIZATION"),
                        selectInput(
                            "p373id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p373plot", "Download")
                    ),
                    mainPanel(plotOutput("p373plot", width = "900px", height = "450px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TOP MARKER EXPRESSIONS IN EACH CELL TYPE"),
                        selectInput(
                            "p38id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p38plot", "Download")
                    ),
                    mainPanel(plotOutput("p38plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
            ),
            tabPanel(
                "4.Sample Pathway Analysis",
                value = "section4",
                h2("Sample Pathway Analysis"),
                h5("Plots showing enrichment analysis and pathway analysis of individual samples."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("ENRICHMENT TERMS OF TOP GENES BAR PLOT"),
                        selectInput(
                            "p26id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p262id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p26plot", "Download")
                    ),
                    mainPanel(plotOutput("p26plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(width = 12,
                                 h3("ENRICHMENT PATHWAYS"),
                                 selectInput(
                                     "p27id",
                                     label = "Choose a sample to display",
                                     choices = "",
                                     selected = ""
                                 )
                    ),
                    mainPanel(fluidRow(
                        column(12, DTOutput('p27table'))
                    ), width = 12),
                    position = "left",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("ENRICHMENT TERMS OF TOP GENES DOT PLOT"),
                        selectInput(
                            "p28id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p282id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p28plot", "Download")
                    ),
                    mainPanel(plotOutput("p28plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("GENE CONCEPT NETWORK CLUSTER"),
                        selectInput(
                            "p29id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p292id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p29plot", "Download")
                    ),
                    mainPanel(plotOutput("p29plot", width = "900px", height = "300px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("NETWORK ENRICHMENT TERMS"),
                        selectInput(
                            "p30id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p302id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p30plot", "Download")
                    ),
                    mainPanel(plotOutput("p30plot", width = "900px", height = "300px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("ENRICHMENT TERMS OF TOP GENES HEATMAP"),
                        selectInput(
                            "p31id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p312id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p31plot", "Download")
                    ),
                    mainPanel(plotOutput("p31plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("ENRICHMENT MAP OF TOP GENES"),
                        selectInput(
                            "p32id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p322id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p32plot", "Download")
                    ),
                    mainPanel(plotOutput("p32plot", width = "900px", height = "1000px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("KEGG ENRICHMENT OF TOP GENES"),
                        selectInput(
                            "p33id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        # sliderInput("p332id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p33plot", "Download")
                    ),
                    mainPanel(plotOutput("p33plot", width = "900px", height = "1000px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("ENRICHMENT TERMS UPSET PLOT"),
                        selectInput(
                            "p34id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p342id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p34plot", "Download")
                    ),
                    mainPanel(plotOutput("p34plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("GSEA GENE ENRICHMENT"),
                        selectInput(
                            "p35id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p352id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p35plot", "Download")
                    ),
                    mainPanel(plotOutput("p35plot", width = "900px", height = "500px")),
                    position = "right",
                )
            ),
            tabPanel(
                "5.Sample Pseudo-time",
                value = "section5",
                h2("Sample Pseudo-time Analysis"),
                h5("Plots showing enrichment analysis and pathway analysis of individual samples."),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE TRAJECTORY ANALYSIS BY CLUSTER"),
                        selectInput(
                            "p73id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                    ),
                    mainPanel(plotlyOutput("p73plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TRAJECTORY ANALYSIS BY CLUSTER - SAMPLE"),
                        selectInput(
                            "p72id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p72plot", "Download")
                    ),
                    mainPanel(plotOutput("p72plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE TRAJECTORY ANALYSIS BY CELL TYPE"),
                        selectInput(
                            "p75id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        )
                    ),
                    mainPanel(plotlyOutput("p75plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TRAJECTORY ANALYSIS BY CELL TYPE - SAMPLE"),
                        selectInput(
                            "p74id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p74plot", "Download")
                    ),
                    mainPanel(plotOutput("p74plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("TRAJECTORY ANALYSIS BY PARTITION - SAMPLE"),
                        selectInput(
                            "p76id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p76plot", "Download")
                    ),
                    mainPanel(plotOutput("p76plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("PSEUDO-TIME ANALYSIS CLUSTER - SAMPLE"),
                        selectInput(
                            "p77id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        sliderInput("p772id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p77plot", "Download")
                    ),
                    mainPanel(plotOutput("p77plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("PSEUDO-TIME ANALYSIS CELL TYPE - SAMPLE"),
                        selectInput(
                            "p79id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        selectInput("p792id",label = "Choose a cell type to display", 
                                    choices = "",
                                    selected = ""),
                        actionButton("p79plot", "Download")
                    ),
                    mainPanel(plotOutput("p79plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("PSEUDO TIME VS PARAMETERS IN CLUSTERS"),
                        selectInput(
                            "p78id",
                            label = "Choose a sample to display",
                            choices = "",
                            selected = ""
                        ),
                        selectInput(
                            "p782id",
                            label = "Choose a parameter to display",
                            choices = "",
                            selected = ""
                        ),
                        selectInput(
                            "p783id",
                            label = "Choose a dimension type to display",
                            choices = "",
                            selected = ""
                        ),
                        actionButton("p78plot", "Download")
                    ),
                    mainPanel(plotOutput("p78plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                # sidebarLayout(
                #     sidebarPanel(
                #         h3("PSEUDO TIME VS PARAMETERS IN CELL TYPES"),
                #         selectInput(
                #             "p84id",
                #             label = "Choose a sample to display",
                #             choices = "",
                #             selected = ""
                #         ),
                #         selectInput(
                #             "p842id",
                #             label = "Choose a parameter to display",
                #             choices = "",
                #             selected = ""
                #         ),
                #         selectInput(
                #             "p843id",
                #             label = "Choose a dimension type to display",
                #             choices = "",
                #             selected = ""
                #         ),
                #         actionButton("p84plot", "Download")
                #     ),
                #     mainPanel(plotOutput("p84plot", width = "900px", height = "600px")),
                #     position = "right",
                # )
            ),
            tabPanel(
                "6.Integrative QC",
                value = "section6",
                h2("Sample Integration Quality Control"),
                h5("Plots showing before and after integration implementation."),
                br(),
                h3("UMAP: BEFORE AND AFTER INTEGRATION - BY SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p39plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p39plot", "Download"),
                br(),
                br(),
                h3("tSNE: BEFORE AND AFTER INTEGRATION - BY SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p40plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p40plot", "Download"),
                br(),
                br(),
                h3("PCA: BEFORE AND AFTER INTEGRATION - BY SAMPLE",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p41plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p42plot", "Download"),
                br(),
                br(),
                h3("UMAP: BEFORE AND AFTER INTEGRATION - BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p42plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p42plot", "Download"),
                br(),
                br(),
                h3("tSNE: BEFORE AND AFTER INTEGRATION - BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p43plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p43plot", "Download"),
                br(),
                br(),
                h3("PCA: BEFORE AND AFTER INTEGRATION - BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p44plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p44plot", "Download"),
                br(),
                br(),
                h3("UMAP: BEFORE AND AFTER INTEGRATION - BY BATCH",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p45plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p45plot", "Download"),
                br(),
                br(),
                h3("tSNE: BEFORE AND AFTER INTEGRATION - BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p46plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p46plot", "Download"),
                br(),
                br(),
                h3("PCA: BEFORE AND AFTER INTEGRATION - BY GROUP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p47plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p47plot", "Download"),
            ),
            tabPanel(
                "7.Integrative Phenotypes",
                value = "section7",
                h2("Integrative Phenotype Discoveries"),
                h5("Plots showing phenotypic discoveries of integrated data."),
                br(),
                h3("INTERACTIVE 3D UMAP CELL TYPE - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotlyOutput("p91plot", width = "1400px", height = "800px"),
                br(),
                br(),
                h3("INTERACTIVE 3D tSNE CELL TYPE - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotlyOutput("p92plot", width = "1400px", height = "800px"),
                br(),
                br(),
                h3("INTERACTIVE 3D UMAP CLUSTER - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotlyOutput("p93plot", width = "1400px", height = "800px"),
                br(),
                br(),
                h3("INTERACTIVE 3D tSNE CLUSTER - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotlyOutput("p94plot", width = "1400px", height = "800px"),
                br(),
                br(),
                br(),
                h3("TOP GENE IN EACH CLUSTER - UMAP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p53plot", width = "1200px", height = "1400px"),
                br(),
                actionButton("p53plot", "Download"),
                br(),
                br(),
                h3("TOP GENE IN EACH CLUSTER - HEATMAP",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p55plot", width = "1200px", height = "1400px"),
                br(),
                actionButton("p55plot", "Download"),
                br(),
                br(),
                h3("CELL TYPE DISCOVERY - UMAP ON INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p56plot", width = "1400px", height = "1200px"),
                br(),
                actionButton("p56plot", "Download"),
                br(),
                br(),
                h3("CELL TYPE DISCOVERY - tSNE ON INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p57plot", width = "1400px", height = "1200px"),
                br(),
                actionButton("p57plot", "Download"),
                br(),
                br(),
                h3("TOP EXPRESSED MARKERS IN CELL TYPES - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p66plot", width = "1400px", height = "1000px"),
                br(),
                actionButton("p66plot", "Download"),
                br(),
                br(),
                h3("TOP EXPRESSED GENES IN CELL TYPES WITH UMAP - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p652plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p652plot", "Download"),
                br(),
                br(),
                h3("TOP EXPRESSED GENES IN CELL TYPES WITH tSNE - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p65plot", width = "1400px", height = "700px"),
                br(),
                actionButton("p65plot", "Download"),
                br(),
                br(),
                h3("MEDIAN EXPRESSION OF TOP MARKERS IN CELL TYPES - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p64plot", width = "1200px", height = "1400px"),
                br(),
                actionButton("p64plot", "Download"),
                br(),
                br(),
                h3("CELL TYPE PROPORTION IN EACH SAMPLE - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p63plot", width = "1400px", height = "1000px"),
                br(),
                actionButton("p63plot", "Download"),
                br(),
                br(),
                h3("SAMPLE PROPORTION IN EACH CELL TYPE - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p67plot", width = "1200px", height = "1200px"),
                br(),
                actionButton("p67plot", "Download"),
                br(),
                br(),
                h3("TOP PC FEATURES - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p48plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p48plot", "Download"),
                br(),
                br(),
                h3("TOP FEATURES IN PC COMPONENTS - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p49plot", width = "1200px", height = "1600px"),
                br(),
                actionButton("p49plot", "Download"),
                br(),
                br(),
                h3("GROUP COMPARISON - UMAP ON INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p50plot", width = "1200px", height = "700"),
                br(),
                actionButton("p50plot", "Download"),
                br(),
                br(),
                h3("GROUP COMPARISON - tSNE ON INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p51plot", width = "1200px", height = "700px"),
                br(),
                actionButton("p51plot", "Download"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("VIOLIN PLOT ON TOP EXPRESSED MARKERS IN INTEGRATED DATA"),
                        sliderInput("p54id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p54plot", "Download")
                    ),
                    mainPanel(plotOutput("p54plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),
                h3("TOP EXPRESSED GENES WITH UMAP - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p602plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p602plot", "Download"),
                br(),
                br(),
                h3("TOP EXPRESSED GENES WITH tSNE - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p60plot", width = "1400px", height = "600px"),
                br(),
                actionButton("p60plot", "Download"),
                br(),
                br(),
                h3("TOP EXPRESSED MARKERS IN CLUSTERS - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p61plot", width = "1400px", height = "1000px"),
                br(),
                actionButton("p61plot", "Download"),
                br(),
                br(),
                h3("SAMPLE PROPORTION IN EACH CLUSTER - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p62plot", width = "1200px", height = "1200px"),
                br(),
                actionButton("p62plot", "Download"),
                br(),
                br(),
                h3("MEDIAN EXPRESSION OF TOP MARKERS IN CLUSTERS - INTEGRATED DATA",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p59plot", width = "1200px", height = "1600px"),
                br(),
                actionButton("p59plot", "Download")
            ),
            tabPanel(
                "8.Integrative PSEUDO-TIME",
                value = "section8",
                h2("Integrative PSEUDO-TIME Analysis"),
                h5("Plots showing PSEUDO-TIME analysis of integrated data."),
                br(),
                h3("TRAJECTORY ANALYSIS BY CELL TYPE - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotlyOutput("p88plot", width = "1400px", height = "900px"),
                br(),
                br(),
                h3("TRAJECTORY ANALYSIS BY CELL TYPE - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p87plot", width = "1400px", height = "900px"),
                br(),
                actionButton("p87plot", "Download"),
                br(),
                br(),
                h3("TRAJECTORY ANALYSIS BY CLUSTER - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotlyOutput("p86plot", width = "1400px", height = "900px"),
                br(),
                br(),
                h3("TRAJECTORY ANALYSIS BY CLUSTER - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p85plot", width = "1400px", height = "900px"),
                br(),
                actionButton("p85plot", "Download"),
                br(),
                br(),
                h3("TRAJECTORY ANALYSIS BY PARTITION - INTEGRATED",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p89plot", width = "1400px", height = "1000px"),
                br(),
                actionButton("p89plot", "Download"),
                br(),
                br(),
                sidebarLayout(
                    sidebarPanel(
                        h3("INTERACTIVE TRAJECTORY ANALYSIS BY CLUSTER"),
                        sliderInput("p90id",label = "Choose a cluster to display", min = 0, max = 10, value = 0, step = 1),
                        actionButton("p90plot", "Download")
                    ),
                    mainPanel(plotOutput("p90plot", width = "900px", height = "600px")),
                    position = "right",
                ),
                br(),
                br(),

            )
        )
    )))
