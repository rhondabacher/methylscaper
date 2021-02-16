ui <- navbarPage("methylscaper",
    tabPanel("Single-cell",
          navbarPage("",
                   tabPanel("Preprocessing",
                            shinyDirButton('folder', 'Select a folder', 'Please select a folder'),
                            verbatimTextOutput("sc_folder_name"), ## we want to display the path
                            textInput("chromosome.number", "Chromosome number"),
                            actionButton("run.subset", label = "Run"),
                            downloadButton("sc_preprocessing_down", label = "Download processed data")),               
                     tabPanel("Seriation",
                              sidebarLayout(
                                sidebarPanel(
                                  fileInput("sc_rds_file", label = "RDS file"),
                                  uiOutput("startPos"),
                                  uiOutput("endPos"),
                                  uiOutput("positionSlider"),
                                  selectInput("sc_ser_method", label = "Seriation Method:",
                                              choices = c("PCA", "ARSA")),
                                  selectInput("sc_refine_method", label = "Refinement Method:",
                                              choices = c("PCA", "HC_average")),
                                  radioButtons("sc_brush_choice", label = "Brushing for:",
                                               choices = c("Refinement", "Weighting"), selected = "Weighting"),
                                  actionButton("sc_force_reverse", label = "Force Reverse"),
                                  verbatimTextOutput("sc_info")
                                ),
                        mainPanel(
                            fluidRow(
                                column(width = 8,
                                    plotOutput(outputId = "sc_seqPlot",brush = "sc_plot_brush",  width = "100%")),
                                column(width = 2, align='left',
                                    selectInput("sc_plot_filetype", label = "Choose file type for saving heatmap", choices = c("PNG", "SVG", "PDF")),
                                    downloadButton("sc_plot_down", label = "Download the heatmap"),
                                    downloadButton("sc_log_down", label = "Download ordering log"))
                                    )
                        )
                        )),
                    tabPanel("Summary Statistics",
                            radioButtons("sc_proportion_choice", label = "Proportion of:", choices = c("Yellow", "Red"), selected = "Yellow"),
                            splitLayout(cellWidths = c("50%", "50%"),
                                plotOutput(outputId = "sc_proportion_color_histogram"),
                                plotOutput(outputId = "sc_percent_C")),
                            splitLayout(cellWidths = c("50%", "50%"),
                              downloadButton("sc_proportion_hist_download", label = "Download histogram"),
                              downloadButton("sc_percentC_plot_download", label = "Download plot")),
                            splitLayout(cellWidths = c("50%", "50%"),
                              downloadButton("sc_proportion_data_download", label = "Download proportion data"),
                              downloadButton("sc_percentC_data_download", label = "Download percentage data"))))),
    tabPanel("Single-molecule",
        navbarPage("",
           tabPanel("Preprocessing",
                    fileInput("fasta.file", label = "FASTA File", accept= c(".fa", ".txt", ".fasta")),
                    fileInput("ref.file", label = "Reference File", accept= c(".fa", ".txt", ".fasta")),
                    #textInput("gch.file.name", label = "GCH File Name"),
                    #textInput("hcg.file.name", label = "HCG File Name"),
                    # textInput("processing.log.name", label = "Processing Log File Name"),
                    actionButton("run.align", label = "Run"),
                    downloadButton("processing_log", label = "Download log file"),
                    downloadButton("sm_preprocessing_down", label = "Download processed data")),
           tabPanel( "Seriation",
                     sidebarLayout(
                       sidebarPanel(
                         fileInput("sm_rds_file", label = "RDS file input", accept= c(".rds", ".RDS")),
                         selectInput("sm_ser_method", label = "Seriation Method:",
                                     choices = c("PCA", "ARSA")),
                         selectInput("sm_refine_method", label = "Refinement Method:",
                                     choices = c("PCA", "HC_average")),
                         radioButtons("sm_brush_choice", label = "Brushing for:",
                                      choices = c("Refinement", "Weighting"), selected = "Weighting"),
                         actionButton("sm_force_reverse", label = "Force Reverse"),
                         verbatimTextOutput("sm_info")
                       ),

                       mainPanel(
                         fluidRow(column(width = 8,
                                         plotOutput(outputId = "sm_seqPlot",
                                                    brush = "sm_plot_brush",  width = "100%")),
                                  column(width = 2, align='left',
                                         selectInput("sm_filetype", label = "File type", choices = c("PNG", "SVG", "PDF")),
                                         downloadButton("sm_plot_down", label = "Download the plot"),
                                         downloadButton("sm_log_down", label = "Download changes log"))

                         )
                       )
                     )),
           tabPanel("Summary Statistics",
                    radioButtons("sm_proportion_choice", label = "Proportion of:", choices = c("Yellow", "Red"), selected = "Yellow"),
                    splitLayout(cellWidths = c("50%", "50%"),
                        plotOutput(outputId = "sm_proportion_color_histogram"),
                        plotOutput(outputId = "sm_percent_C")),
                    splitLayout(cellWidths = c("50%", "50%"),
                      downloadButton("sm_proportion_hist_download", label = "Download histogram"),
                      downloadButton("sm_percentC_plot_download", label = "Download plot")),
                    splitLayout(cellWidths = c("50%", "50%"),
                      downloadButton("sm_proportion_data_download", label = "Download proportion data"),
                      downloadButton("sm_percentC_data_download", label = "Download percentage data"))))))

