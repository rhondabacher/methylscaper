ui <- navbarPage("methylScaper",
                 tabPanel("Single-cell",
                          navbarPage("",
                                     tabPanel("Seriation",
                                              sidebarLayout(
                                                sidebarPanel(
                                                  fileInput("gch_seq_file", label = "GCH Sequence RDS file"),
                                                  fileInput("hcg_seq_file", label = "HCG Sequence RDS file"),
                                                  numericInput("startPos", label = "Starting position",value = 105636488, min = 1),
                                                  numericInput("endPos", label = "Ending position", value = 105636993, min = 1),
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
                                                  fluidRow(column(width = 8,
                                                                  plotOutput(outputId = "sc_seqPlot",
                                                                             brush = "sc_plot_brush",  width = "100%")),
                                                           column(width = 2, align='left',
                                                                  selectInput("sc_plot_filetype", label = "File type", choices = c("PNG", "SVG", "PDF")),
                                                                  downloadButton("sc_plot_down", label = "Download the plot"),
                                                                  downloadButton("sc_log_down", label = "Download changes log"))

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
                                              fileInput("fasta.file", label = "FASTA File"),
                                              fileInput("ref.file", label = "Reference File"),
                                              textInput("gch.file.name", label = "GCH File Name"),
                                              textInput("hcg.file.name", label = "HCG File Name"),
                                              textInput("processing.log.name", label = "Processing Log File Name"),
                                              actionButton("run.align", label = "Run")),
                                     tabPanel( "Seriation",
                                               sidebarLayout(
                                                 sidebarPanel(
                                                   fileInput("sm_gch_file", label = "GCH Data file input"),
                                                   fileInput("sm_hcg_file", label = "HCG Data file input"),
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
