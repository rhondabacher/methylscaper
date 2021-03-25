library(shinyFiles)
library(shiny)
library(shinyjs)

ui <- navbarPage("methylscaper",

    tabPanel("Single-cell",
          navbarPage("",
                   tabPanel("Preprocessing",
                   
                   fluidRow(
                       column(3,
                            tags$br(),
                            tags$br(),
                            shinyDirButton('folder', 'Select a folder', 'Please select a folder'),
                            verbatimTextOutput("sc_folder_name"), ## display the path
                            textInput("chromosome_number", "Chromosome number"),
                            actionButton("run_subset", label = "Run"),
                            ),
                        column(6, offset = 1,
                            tags$h3("Directions"),
                            tags$ul(
                                tags$li("Select the folder containing the summarized data from the single-cell experiment. This summarized data
                                should be in the form of three columns (chromosome, position, methylation status). This is the typical
                                format obtained from programs such as Bismarck. The folder should contain two subfolders with 
                                the names 'met' and 'acc' which contain the files with the endogenous methylation and the 
                                accessibility methylation, respectively."),
                                tags$li("As these files are large, we preprocess a single chromosome at a time. 
                                    Indicate the desired chromosome as a number or as X, Y, MT."),
                                tags$li("After pressing 'Run', a progress bar will appear in the bottom right.")),
                                tags$br(),
                                tags$p("After processing and subsetting the data for visualization, you may download the data below. 
                                The data will be downloaded in RDS format and can be loaded into the Seriation tab."),   
                            
                            downloadButton("sc_preprocessing_down", label = "Download Processed Data")
                            )),
                    tags$br(),
                    tags$br(),
                    fluidRow(
                        tags$p("For additional information or to view Frequently Asked Questions, see the methylscaper 
                        vignette:", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")),
        
                        tags$p("Questions can also be submitted on our GitHub page:", 
                        tags$a(href="https://github.com/rhondabacher/methylscaper/issues", "Report issue/bug")))
                   ),        
                             
                    
                tabPanel("Seriation",
                              
                    fluidRow(
                         column(11, 
                         tags$p("Upload the RDS file obtained in the Preprocessing tab below in the input labelled 
                         'RDS File Input'."),
                         tags$p("To move along the genome, we have pre-loaded gene locations for Human and Mouse 
                         for the chromosome selected in the Preprocessing tab. Select a gene and then a slider will appear 
                         to refine the genomic location."),
                         tags$p("With Brushing set to Weighting, drag the mouse/cursor to select basepairs (columns) 
                         on which to weight the ordering algorithm. Once the weighting is set, change Brushing to 
                         refinement and highlight the cells (rows) to refine the ordering."),
                         tags$p("Detailed
                         descriptions of the various options for seriation, weighting, and refinement are in the methylscaper 
                         vignette: ", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")))
                     ),
                     fluidRow(sidebarLayout(
                                sidebarPanel(
                                  actionButton("sc_demo_data", label = "Load Example Data"),
                                  fileInput("sc_rds_file", label = "RDS File Input"),
                                  radioButtons("organism_choice", label = "Choose Organism:", choices = c("Human", "Mouse", "Other"), 
                                                          selected = character(0)),
                                  selectizeInput("geneList", label="Choose Gene", choices=NULL,
                                  options = list(
                                      placeholder = 'Begin typing or select from below',
                                      onInitialize = I('function() { this.setValue(""); }')
                                    )),                        
                                  uiOutput("startPos"),
                                  uiOutput("endPos"),
                                  uiOutput("positionSlider"),
                                  selectInput("sc_ser_method", label = "Seriation Method:",
                                              choices = c("PCA", "ARSA")),
                                  selectInput("sc_refine_method", label = "Refinement Method:",
                                              choices = c("PCA", "HC_average")),
                                  radioButtons("sc_brush_choice", label = "Brushing for:",
                                               choices = c("Refinement", "Weighting"), selected = "Weighting"),
                                  actionButton("sc_force_reverse", label = "Reverse Current Ordering"),
                                  verbatimTextOutput("sc_info"),
																	width = 3
                                ),
                        mainPanel(
                            fluidRow(
                                column(width = 10,
																	  useShinyjs(),
                                    plotOutput(outputId = "sc_seqPlot",brush = "sc_plot_brush",  width = "100%")),
                                column(width = 2, align='left',
                                    selectInput("sc_plot_filetype", label = "Choose file type for saving heatmap", choices = c("PNG", "PDF")),
                                    downloadButton("sc_plot_down", label = "Download Heatmap"),
                                    downloadButton("sc_log_down", label = "Download Ordering Log"))
                                    )
                        )
                        ))),
                    tabPanel("Summary Statistics",
                    fluidRow(
                        column(10, 
                        tags$p("Below are summary plots of the methylation data. The histogram on the left
                        calulates the proportion of methylated basepairs in each cell. Choose whether to calculate this for
                        the endogenous methylation or the accessibility methylation. The plot on the right calculates the percent of
                        cells having methylation along the region, and the dots indicate either GCH (yellow) or HCG (red) sites."
                        )
                       )),
                      fluidRow( radioButtons("sc_proportion_choice", label = "Proportion of:", 
                                      choices = c("Accessibility Methylation", "Endogenous Methylation"), selected = "Accessibility Methylation"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                plotOutput(outputId = "sc_proportion_color_histogram"),
                                plotOutput(outputId = "sc_percent_C")),
                            splitLayout(cellWidths = c("50%", "50%"),
                              downloadButton("sc_proportion_hist_download", label = "Download Histogram"),
                              downloadButton("sc_percentC_plot_download", label = "Download Plot")),
                            splitLayout(cellWidths = c("50%", "50%"),
                              downloadButton("sc_proportion_data_download", label = "Download Proportion Data"),
                              downloadButton("sc_percentC_data_download", label = "Download Percentage Data")))))),
    tabPanel("Single-molecule",
        navbarPage("",
           tabPanel("Preprocessing",
           fluidRow(
               column(3,
                    tags$br(),
                    tags$br(),
                    fileInput("fasta_file", label = "FASTA File", accept= c(".fa", ".txt", ".fasta")),
                    fileInput("ref_file", label = "Reference File", accept= c(".fa", ".txt", ".fasta")),
                    actionButton("run_align", label = "Run"),
                    ),
                column(6, offset = 1,
                    tags$h3("Directions"),
                    tags$ul(
                        tags$li("Upload the reads or sequences from the single-molecule experiment in FASTA format in the input labelled
                                'FASTA File'."),
                        tags$li("Upload the reference sequences for the gene or genomic location of interest in the input 
                                labelled 'Reference File', this should also be in FASTA format."),
                        tags$li("After pressing 'Run', a progress bar will appear in the bottom right.")),
                        tags$br(),
                        tags$p("After aligining and processing the reads for visualization, you may download the data and the 
                    preprocessing log below. The data will be downloaded in RDS format and can be loaded into the Seriation tab. 
                    The processing log contains details on the read alignments."),
                    downloadButton("sm_preprocessing_down", label = "Download Processed Data"),
                    downloadButton("processing_log", label = "Download Log File")
            )),
            tags$br(),
            tags$br(),
            fluidRow(
                tags$p("For additional information or to view Frequently Asked Questions, see the methylscaper 
                vignette:", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")),
                
                tags$p("Questions can also be submitted on our GitHub page:", 
                tags$a(href="https://github.com/rhondabacher/methylscaper/issues", "Report issue/bug")))
           ),
           tabPanel( "Seriation",
                     fluidRow(
                         column(11, 
                         tags$p("Upload the RDS file obtained in the Preprocessing tab below in the input labelled 
                         'RDS File Input'."),
                         tags$p("With Brushing set to Weighting, drag the mouse/cursor to select basepairs (columns) 
                         on which to weight the ordering algorithm. Once the weighting is set, change Brushing to 
                         refinement and highlight the molecules (rows) to refine the ordering."),
                         tags$p("Detailed
                         descriptions of the various options for seriation, weighting, and refinement are in the methylscaper 
                         vignette: ", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")))
                     ),
                     fluidRow(sidebarLayout(
                       sidebarPanel(
                         actionButton("sm_demo_data", label = "Load Example Data"),
                         fileInput("sm_rds_file", label = "RDS File Input", accept= c(".rds", ".RDS")),
                         selectInput("sm_ser_method", label = "Seriation Method:",
                                     choices = c("PCA", "ARSA")),
                         selectInput("sm_refine_method", label = "Refinement Method:",
                                     choices = c("PCA", "HC_average")),
                         radioButtons("sm_brush_choice", label = "Brushing for:",
                                      choices = c("Refinement", "Weighting"), selected = "Weighting"),
                         actionButton("sm_force_reverse", label = "Reverse Current Ordering"),
                         verbatimTextOutput("sm_info"),
												 width = 3
                       ),

                       mainPanel(
                         fluidRow(column(width = 10,
                                         plotOutput(outputId = "sm_seqPlot",
                                                    brush = "sm_plot_brush",  width = "100%")),
                                  column(width = 2, align='left',
                                         selectInput("sm_filetype", label = "File type", choices = c("PNG", "PDF")),
                                         downloadButton("sm_plot_down", label = "Download Heatmap"),
                                         downloadButton("sm_log_down", label = "Download Ordering Log"))
                         )
                       )
                     ))),
           tabPanel("Summary Statistics",
               fluidRow(
                   column(10, 
                   tags$p("Below are summary plots of the methylation data. The histogram on the left
                   calulates the proportion of methylated basepairs in each molecule. Choose whether to calculate this for
                   the endogenous methylation or the accessibility methylation. The plot on the right calculates the percent of
                   molecules having methylation along the region, and the dots indicate either GCH (yellow) or HCG (red) sites."
                   )
                  )),
                 fluidRow( radioButtons("sm_proportion_choice", label = "Proportion of:", 
                                 choices = c("Accessibility Methylation", "Endogenous Methylation"), selected = "Accessibility Methylation"),
                    splitLayout(cellWidths = c("50%", "50%"),
                        plotOutput(outputId = "sm_proportion_color_histogram"),
                        plotOutput(outputId = "sm_percent_C")),
                    splitLayout(cellWidths = c("50%", "50%"),
                      downloadButton("sm_proportion_hist_download", label = "Download Histogram"),
                      downloadButton("sm_percentC_plot_download", label = "Download Plot")),
                    splitLayout(cellWidths = c("50%", "50%"),
                      downloadButton("sm_proportion_data_download", label = "Download Proportion Data"),
                      downloadButton("sm_percentC_data_download", label = "Download Percentage Data")))))))
