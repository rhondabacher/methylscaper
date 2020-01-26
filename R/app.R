
ui <- navbarPage("methylScaper",
                 tabPanel("Preprocessing",
                          fileInput("fasta.file", label = "FASTA File"),
                          fileInput("ref.file", label = "Reference File"),
                          textInput("gch.file.name", label = "GCH File Name"),
                          textInput("hcg.file.name", label = "HCG File Name"),
                          actionButton("run.align", label = "Run")),
                 tabPanel("Sequence Plot",

    sidebarLayout(
        sidebarPanel(
            fileInput("gch.file", label = "GCH Data file input (csv)"),
            fileInput("hcg.file", label = "HCG Data file input (csv)"),
            selectInput("method", label = "Seriation Method:", 
                        choices = c("PCA", "ARSA")),
            selectInput("refineMethod", label = "Refinement Method:", 
                        choices = c("PCA", "HC_average")),
            radioButtons("brush.choice", label = "Brushing for:",
                               choices = c("Refinement", "Weighting"), selected = "Weighting"),
            actionButton("force.reverse", label = "Force Reverse"),
             verbatimTextOutput("info")
        ),

        mainPanel(
         fluidRow(column(width = 8,
                plotOutput(outputId = "seqPlot",
                           brush = "plot_brush",  width = "100%")),
                 column(width = 2, align='left',       
                selectInput("filetype", label = "File type", choices = c("PNG", "SVG", "PDF")),
                downloadButton("down", label = "Download the plot"),
                downloadButton("down_log", label = "Download changes log"))
          
        )
       )
  ))
)

server <- function(input, output) {
  
  
    # alignment handling
    observeEvent(input$run.align, { 
        ref <- read.fasta(input$ref.file$datapath)
        fasta <- read.fasta(input$fasta.file$datapath)
        
        progress <- Progress$new()
        progress$set(message = "Beginning alignment", value = 0)
        on.exit(progress$close())
        
        updateProgress <- function(value = NULL, message = NULL, detail = NULL) { 
          progress$set(value = value, message = message, detail = detail)}
        
        align.out <- runAlign(ref, fasta, fasta.subset = (1:100), updateProgress = updateProgress)
        write.table(align.out$hcg, file=input$hcg.file.name, quote=F, row.names = F, sep="\t")
        write.table(align.out$gch, file=input$gch.file.name, quote=F, row.names = F, sep="\t")
    })
  
  
  
  
    
    actionsLog <- reactiveValues(log = c("Loading Day7 data"))
    input.Data <- reactiveValues(gch = methylScaper::day7$gch, hcg = methylScaper::day7$hcg)
    
    observe({if (!is.null(input$gch.file) & !is.null(input$hcg.file))
    {
        temp.gch <- read.table(input$gch.file$datapath, header=T, row.names = 1, 
                                 stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
        temp.hcg <- read.table(input$hcg.file$datapath, header=T, row.names = 1, 
                                     stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
        if (nrow(temp.gch) == nrow(temp.hcg))
        {
          coordinatesObject$refine.start <- 0
          coordinatesObject$refine.stop <- 0
          coordinatesObject$weight.start <- 0
          coordinatesObject$weight.stop <- 0
          input.Data$gch <- temp.gch
          input.Data$hcg <- temp.hcg
          isolate({
            actionsLog$log <- c(actionsLog$log, paste("Loading GCH file:", input$gch.file$name))
            actionsLog$log <- c(actionsLog$log, paste("Loading HCG file:", input$hcg.file$name))
          })
        }
        
    }})
    
    # this object keeps track of the coordinates for refinement and weighting
    coordinatesObject <- reactiveValues(refine.start = 0, refine.stop = 0, 
                                        weight.start = 0, weight.stop = 0, weight.color = "red")
    # now construct the orderObject
    orderObject <- reactiveValues(toClust = 0, order1 = 0)
    observe({
        tempObj <- buildOrderObjectShiny(input.Data$gch, input.Data$hcg, input$method, coordinatesObject)
        orderObject$order1 <- tempObj$order1
        orderObject$toClust <- tempObj$toClust
        isolate({
            actionsLog$log <- c(actionsLog$log, 
                                paste("Ordering with", input$method))
            })
        })
    
    # this handles updates to coordinatesObject
    observeEvent(input$plot_brush, {
            n <- nrow(input.Data$gch)
            m <- ncol(input.Data$hcg)
            processed.brush <- handleBrushCoordinates(input$plot_brush, n, m)
            
            if (isolate(input$brush.choice) == "Weighting")
            {
                coordinatesObject$refine.start <- 0
                coordinatesObject$refine.stop <- 0
                coordinatesObject$weight.start <- processed.brush$first.col
                coordinatesObject$weight.stop <- processed.brush$last.col
                coordinatesObject$weight.color <- processed.brush$weight.color
                isolate({
                    actionsLog$log <- c(actionsLog$log, 
                                        paste("Weighting", processed.brush$weight.color, "columns", 
                                               processed.brush$first.col, "to",
                                              processed.brush$last.col))
                })
            }
            if (isolate(input$brush.choice) == "Refinement")
            {
                s <- processed.brush$first.row
                f <- processed.brush$last.row
                if (s != f)
                {
                  coordinatesObject$refine.start <- s
                  coordinatesObject$refine.stop <- f
                    orderObject$order1 <- refineOrderShiny(isolate(orderObject), 
                                                           refine.method = isolate(input$refineMethod), 
                                                           coordinatesObject)
                    isolate({
                        actionsLog$log <- c(actionsLog$log, 
                                           paste("Refining rows", 
                                                 processed.brush$first.row, "to",
                                                 processed.brush$last.row))
                        actionsLog$log <- c(actionsLog$log, 
                                            paste("Applying refinement with", input$refineMethod)) 
                    })
                }
            }
            
            
    })
    
    observeEvent( input$force.reverse, {
      isolate({
        if (coordinatesObject$refine.start == coordinatesObject$refine.stop) 
          {
            orderObject$order1 <- rev(orderObject$order1)
            actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to", nrow(input.Data$gch)))
          }
        else 
          {
              orderObject$order1[coordinatesObject$refine.start : coordinatesObject$refine.stop] <-
              orderObject$order1[coordinatesObject$refine.stop : coordinatesObject$refine.start]
              actionsLog$log <- c(actionsLog$log, 
                                  paste("Reversing rows", coordinatesObject$refine.start, 
                                        "to", coordinatesObject$refine.stop))
              
        }
      })
    })
    

    
    output$seqPlot <- renderPlot({ 
        obj <- orderObject
        # print(str(as.list(isolate(coordinatesObject))))
        makePlot(obj,isolate(coordinatesObject))
        }, height=600, width=600)
    
    output$down <- downloadHandler(
        filename = function(){
          if (input$filetype == "PNG") return("plot.png")
          if (input$filetype == "SVG") return("plot.svg")
          if (input$filetype == "PDF") return("plot.pdf")
        },
        content = function(file){
            if (input$filetype == "PNG") png(file, res=300)
            if (input$filetype == "SVG") svg(file)
            if (input$filetype == "PDF") pdf(file)
          
            makePlot(orderObject, coordinatesObject, drawLines = FALSE, plotFAST = FALSE)
            dev.off()
        }
    )
    
    output$down_log <- downloadHandler(
        filename = function(){
            "changes.log"
        },
        content = function(file){
            fileConn <- file(file)
            writeLines(actionsLog$log, fileConn)
            close(fileConn)
        }
    )
    
    output$info <- renderText({
        paste0("Refinement selection: ", coordinatesObject$refine.start, " ", coordinatesObject$refine.stop, "\n",
               "Weighting selection: ", coordinatesObject$weight.start, " ", coordinatesObject$weight.stop)
    })
}





# Run the application 
#' methylScaper
#' 
#' Runs the methylScaper Shiny app.
#' 
#' @import shiny
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom utils read.table write.table
#' @export
methylScaper <- function() {
	options(shiny.maxRequestSize = 10000*1024^2) 
	shinyApp(ui = ui, server = server)}
