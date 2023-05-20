library(shiny)
library(shinythemes)
library(ggplot2)
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
          }, 15000000000000)
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
    
    title = "scCNV | SCA",
    theme = shinytheme("cosmo"),
    titlePanel(# strong("Single Cell RNA-Sequencing Online Analysis Platform",align = "center")),
        tags$link(rel="icon", href="DB/favicon"),
        h2(
            HTML("<b>Single Cell CNV Online Analysis Platform</b>"),
            style = "text-align:center"
        )),
    h2(
        "Single Cell CNV Online Analysis Platform",
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
                        'Upload a zip file containing 10X format CNV mappable regions .bed file, cnv calls .bed file and summary metrics .csv file for all your samples (please click submit only after all your files have been uploaded):'),
                    multiple = TRUE,
                    accept = c(
                        '.tar',
                        '.tar.gz',
                        'application/zip',
                        '.zip',
                        '.gz'
                    )
                ),
                downloadButton("downloadExample1", label = "scCNV Input Example"),
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
                downloadButton("downloadExample2", label = "scCNV Metadata Example"),
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
                "1.Ploidy Info",
                value = "section1",
                h2("Ploidy Information"),
                h5("Metrics on single-cell level Copy Number Variations (CNV) of samples. Please waiting patiently while plotting is done."),
                sidebarLayout(
                    sidebarPanel(
                        h3("PLOIDY INFO"),
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
                
                h3("MEAN SINGLE-CELL PLOIDY (SAMPLES)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p2plot", width = "1200px", height = "800px"),
                actionButton("p2plot", "Screenshot"),
                br(),
                br()
            ),
            tabPanel(
                "2.Dimension Reduction and Clustering",
                value = "section2",
                h2("Dimension Reduction and Clustering"),
                h5("Clustering samples together based on their CNV information. Please waiting patiently while plotting is done."),
                br(),
                h3("DAPC COMPONENTS OF CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p3plot", width = "1200px", height = "800px"),
                br(),
                actionButton("p3plot", "Screenshot"),
                br(),
                br(),
                h3("UMAP - CLUSTERS (ALL SAMPLES)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p4plot", width = "1200px", height = "800px"),
                actionButton("p4plot", "Screenshot"),
                br(),
                br(),
                h3("SINGLE-CELL CNV EVENTS BY CHROMOSOME POSITIONS (CLUSTERS)",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p5plot", width = "1200px", height = "800px"),
                actionButton("p5plot", "Screenshot"),
                br(),
                br(),
                h3("PLOIDY INFORMATION BASED ON TOP 50 CNV EVENTS IN DAPC COMPONENTS 1 AND 2",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p6plot", width = "1400px", height = "800px"),
                actionButton("p6plot", "Screenshot"),
                br(),
                br(),
                h3("CLUSTER PROPORTION IN SAMPLES",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p7plot", width = "1400px", height = "800px"),
                actionButton("p7plot", "Screenshot"),
                br(),
                br(),
                h3("PHYLOGENETIC TREE OF COPY NUMBER CLUSTERS",
                   style={'background-color:#EDF0F1;
                       text-align: center; padding: 25px 0;border: 2px solid #D6D9DC;'}),
                plotOutput("p8plot", width = "1400px", height = "1200px"),
                actionButton("p8plot", "Screenshot"),
                br(),
                br()
            )
        )
    )
))
