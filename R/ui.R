library(shinyFiles)
library(shiny)
library(shinyjs)

ui <- navbarPage("methylscaper",id="big_tab",

    tabPanel("Single-cell",
          navbarPage("", id="seriate_sc",
                   tabPanel("Preprocessing",
                   
                   fluidRow(
                       column(4,
                            tags$br(),
                            tags$br(),
														fileInput("sc_met_files", label = "Endogenous Methylation Files", multiple=TRUE),
														fileInput("sc_acc_files", label = "Accessibility Methylation Files", multiple=TRUE),
                            textInput("chromosome_number", "Chromosome number"),
                            actionButton("run_subset", label = "Run"),
                            ),
                        column(6, offset = 1,
                            tags$h3("Directions"),
                            tags$ul(
                                tags$li("Select all methylation files of each type when uploading data from the single-cell experiment. This summarized data
                                should be from Bismark's methylation extractor or in the form of three columns (chromosome, position, methylation rate)."),
                                tags$li("As these files are large, we process a single chromosome at a time. 
                                    Indicate the desired chromosome as a number (e.g., 1 - 22) or as X, Y, MT."),
                                tags$li("After pressing 'Run', a progress bar will appear in the bottom right.")),
                                tags$br(),
                                tags$p("After processing and subsetting the data for visualization, you may download the data below. 
                                The data will be downloaded in RDS format and can be loaded into the Visualization tab."),   
                            
                            shinyjs::disabled(downloadButton("sc_preprocessing_down", label = "Download Processed Data"))
                            )),
                    tags$br(),
                    tags$br(),
                    fluidRow(
                        tags$p("For additional information or to view Frequently Asked Questions, the linked methylscaper 
                vignette is:", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")),
        
                        tags$p("Questions can also be submitted on our GitHub page:", 
                        tags$a(href="https://github.com/rhondabacher/methylscaper/issues", "Report issue/bug")))
                   ),        
                             
                    
                tabPanel("Visualization", 
                              
                    fluidRow(
                         column(11, 
                         tags$p("Upload the RDS file obtained in the Preprocessing tab below in the input labelled 
                         'RDS File Input'. Alternatively, explore the app by selecting 'Load Example Data' below. The example data
												 is from Mouse on Chromosome 19 subset to a 40 kbp region around 8,967,041. Some genes that are located in this region
												 are", tags$em("Eef1g,"), tags$em("Mta2,"), "and", tags$em("Tut1.")),
                         tags$p("To move along the genome, we have pre-loaded gene locations for Human (GRCh38) and Mouse (GRCm39)
                         for the chromosome selected in the Preprocessing tab. Select a gene and then a slider will appear 
                         to refine the genomic location."),
                         tags$p("With Brushing set to Weighting, drag the mouse/cursor to select basepairs (columns) 
                         on which to weight the ordering algorithm. Once the weighting is set, change Brushing to 
                         refinement and highlight the cells (rows) to refine the ordering."),
                         tags$p("Detailed descriptions of the various options for seriation, weighting, and refinement are in the methylscaper 
                         vignette: ", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")))
                     ),
                     fluidRow(sidebarLayout(
                                sidebarPanel(
                                  actionButton("sc_demo_data", label = "Load Example Data", style = "background-color: rgb(209, 238, 238)"),
                                  tags$br(),
                                  tags$br(),
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
                                    selectInput("sc_plot_filetype", label = "Choose file type for saving heatmap", choices = c("PNG", "PDF", "SVG")),
                                    shinyjs::disabled(downloadButton("sc_plot_down", label = "Download Heatmap")),
                                    shinyjs::disabled(downloadButton("sc_log_down", label = "Download Ordering Log")))
                                    )
                        )
                        ))),
		            tabPanel("Summary Statistics",
		                fluidRow(
		                    column(10, 
		                    tags$p("Below are summary plots of the methylation data. The histogram on the left
		                    calulates the proportion of methylated bases within each cell. Choose whether to calculate this for
		                    the endogenous methylation or the accessibility methylation. The plot in the middle calculates the percent of
		                    cells having methylation at each GCH (yellow) or HCG (red) site in the region. The plot on the right calculates
		 									 a population averaged methylation status over the entire region. The window size can be adjusted. This plot
											 is similar to the one in the middle but is easier to interpret when the region is large."
		                    )
		                 )),
		                  fluidRow(column(10, radioButtons("sc_proportion_choice", label = "Proportion of:", 
		                                  choices = c("Accessibility Methylation", "Endogenous Methylation"), selected = "Accessibility Methylation")),
		 													column(2, numericInput("sc_window_choice", label = "Window size:", 
		 													                                  value=20))),
		 						     fluidRow(																								
		                     splitLayout(cellWidths = c("33%", "33%", "33%"),
		                         plotOutput(outputId = "sc_proportion_color_histogram"),
		                         plotOutput(outputId = "sc_percent_C"),
		 												plotOutput(outputId = "sc_avg_c")),
		                     splitLayout(cellWidths = c("33%", "33%", "33%"),
		                       shinyjs::disabled(downloadButton("sc_proportion_hist_download", label = "Download Histogram")),
		                       shinyjs::disabled(downloadButton("sc_percentC_plot_download", label = "Download Plot")),
		 											shinyjs::disabled(downloadButton("sc_avg_c_plot_download", label = "Download Plot"))),
		                     splitLayout(cellWidths = c("33%", "33%", "33%"),
		                       shinyjs::disabled(downloadButton("sc_proportion_data_download", label = "Download Proportion Data")),
		                       shinyjs::disabled(downloadButton("sc_percentC_data_download", label = "Download Percentage Data")),
		 											shinyjs::disabled(downloadButton("sc_avg_c_data_download", label = "Download Averaged Percent Data"))))
            
		 											))),
    tabPanel("Single-molecule",
        navbarPage("",  id="seriate_sm",
           tabPanel("Preprocessing",
           fluidRow(
               column(3,
                    tags$br(),
                    tags$br(),
                    fileInput("fasta_file", label = "Reads .FASTA File", accept= c(".fa", ".txt", ".fasta"), multiple=TRUE),
                    fileInput("ref_file", label = "Reference .FASTA File", accept= c(".fa", ".txt", ".fasta"), multiple=TRUE),
                    actionButton("run_align", label = "Run"),
                    ),
                column(6, offset = 1,
                    tags$h3("Directions"),
                    tags$ul(
                        tags$li("Upload the reads or sequences from the single-molecule experiment in FASTA format in the input labelled
                                'FASTA File'. There is a 10MB upload limit for files on the webserver. For larger files, 
																please run methylscaper locally. Additional details are provided in the vignette linked at the bottom of this page."),
                        tags$li("Upload the reference sequence for the gene or genomic location of interest in the input 
                                labelled 'Reference File', this should also be in FASTA format."),
                        tags$li("After pressing 'Run', a progress bar will appear in the bottom right.")),
                        tags$br(),
                        tags$p("After aligining and processing the reads for visualization, you may download the data and the 
                    preprocessing log below. The data will be downloaded in RDS format and can be loaded into the Visualization tab. 
                    The processing log contains details on the read alignments."),
                    shinyjs::disabled(downloadButton("sm_preprocessing_down", label = "Download Processed Data")),
                    shinyjs::disabled(downloadButton("processing_log", label = "Download Log File"))
            )),
            tags$br(),
            tags$br(),
            fluidRow(
                tags$p("For additional information or to view Frequently Asked Questions, the linked methylscaper 
                vignette is:", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")),
                
                tags$p("Questions can also be submitted on our GitHub page:", 
                tags$a(href="https://github.com/rhondabacher/methylscaper/issues", "Report issue/bug")))
           ),
           tabPanel( "Visualization",
                     fluidRow(
                         column(11, 
                         tags$p("Upload the RDS file obtained in the Preprocessing tab below in the input labelled 
                         'RDS File Input'. Alternatively, explore the app by selecting 'Load Example Data'."),
                         tags$p("With Brushing set to Weighting, drag the mouse/cursor to select basepairs (columns) 
                         on which to weight the ordering algorithm. Once the weighting is set, change Brushing to 
                         refinement and highlight the molecules (rows) to refine the ordering."),
                         tags$p("Detailed
                         descriptions of the various options for seriation, weighting, and refinement are in the methylscaper 
                         vignette: ", tags$a(href="http://www.methylscaper.com/content/vignette.pdf", "PDF")))
                     ),
                     fluidRow(sidebarLayout(
                       sidebarPanel(
                         actionButton("sm_demo_data", label = "Load Example Data", style = "background-color: rgb(209, 238, 238)"),
                         tags$br(),
                         tags$br(),
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
                                         selectInput("sm_filetype", label = "File type", choices = c("PNG", "PDF", "SVG")),
                                         shinyjs::disabled(downloadButton("sm_plot_down", label = "Download Heatmap")),
                                         shinyjs::disabled(downloadButton("sm_log_down", label = "Download Ordering Log")))
                         )
                       )
                     ))),
           tabPanel("Summary Statistics",
               fluidRow(
                   column(10, 
                   tags$p("Below are summary plots of the methylation data. The histogram on the left
                   calulates the proportion of methylated bases within each molecule. Choose whether to calculate this for
                   the endogenous methylation or the accessibility methylation. The plot in the middle calculates the percent of
                   molecules having methylation at each GCH (yellow) or HCG (red) site in the region. The plot on the right calculates
									 a population averaged methylation status over the entire region. The window size can be adjusted. This plot
											 is similar to the one in the middle but is easier to interpret when the region is large."
                   )
                  )),
                 fluidRow(column(10, radioButtons("sm_proportion_choice", label = "Proportion of:", 
                                 choices = c("Accessibility Methylation", "Endogenous Methylation"), selected = "Accessibility Methylation")),
													column(2, numericInput("sm_window_choice", label = "Window size:", 
													                                  value=20))),
						     fluidRow(																								
                    splitLayout(cellWidths = c("33%", "33%", "33%"),
                        plotOutput(outputId = "sm_proportion_color_histogram"),
                        plotOutput(outputId = "sm_percent_C"),
												plotOutput(outputId = "sm_avg_c")),
                    splitLayout(cellWidths = c("33%", "33%", "33%"),
                      shinyjs::disabled(downloadButton("sm_proportion_hist_download", label = "Download Histogram")),
                      shinyjs::disabled(downloadButton("sm_percentC_plot_download", label = "Download Plot")),
											shinyjs::disabled(downloadButton("sm_avg_c_plot_download", label = "Download Plot"))),
                    splitLayout(cellWidths = c("33%", "33%", "33%"),
                      shinyjs::disabled(downloadButton("sm_proportion_data_download", label = "Download Proportion Data")),
                      shinyjs::disabled(downloadButton("sm_percentC_data_download", label = "Download Percentage Data")),
											shinyjs::disabled(downloadButton("sm_avg_c_data_download", label = "Download Averaged Percent Data"))))
                    
											))))
