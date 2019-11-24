options(shiny.maxRequestSize = 10000*1024^2) 
ui <- fluidPage(

    # Application title
    titlePanel("methylScaper"),

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
             verbatimTextOutput("info")
        ),

        mainPanel(
           plotOutput(outputId = "seqPlot",
                      brush = "plot_brush"),
           selectInput("filetype", label = "File type", choices = c("PNG", "SVG")),
           downloadButton("down", label = "Download the plot"),
           downloadButton("down_log", label = "Download changes log")
          
        )
    )
)

server <- function(input, output) {
    
    actionsLog <- reactiveValues(log = c("Loading Day7 data"))
    input.Data <- reactiveValues(gch = day7$gch, hcg = day7$hcg)
    
    observe({if (!is.null(input$gch.file) & !is.null(input$hcg.file))
    {
        input.Data$gch <- read.table(input$gch.file$datapath, header=T, row.names = 1, 
                                 stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
        input.Data$hcg <- read.table(input$hcg.file$datapath, header=T, row.names = 1, 
                                     stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
        isolate({
            actionsLog$log <- c(actionsLog$log, paste("Loading GCH file:", input$gch.file$name))
            actionsLog$log <- c(actionsLog$log, paste("Loading HCG file:", input$hcg.file$name))
        })
    }})
    
    # this object keeps track of the coordinates for refinement and weighting
    coordinatesObject <- reactiveValues(refine.start = 0, refine.stop = 0, 
                                        weight.start = 0, weight.stop = 0, weight.color = "red")
    # now construct the orderObject
    orderObject <- reactiveValues(toClust = 0, order1 = 0)
    observe({
        print("construction orderObjct")
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
            print("updating coordinatesObject")
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
                coordinatesObject$refine.start <- processed.brush$first.row
                coordinatesObject$refine.stop <- processed.brush$last.row
                isolate({
                    actionsLog$log <- c(actionsLog$log, 
                                        paste("Refining rows", 
                                              processed.brush$first.row, "to",
                                              processed.brush$last.row))
                }) 
            }
            
            print("starting to handle refinement")
            s <- coordinatesObject$refine.start
            f <-coordinatesObject$refine.stop
            if (s != 0 & f != 0)
            {
                orderObject$order1 <- refineOrderShiny(isolate(orderObject), 
                                                       refine.method = isolate(input$refineMethod), 
                                                       coordinatesObject)
                print("orderObject updated with refinement")
                isolate({
                    actionsLog$log <- c(actionsLog$log, 
                                        paste("Applying refinement with", input$refineMethod))
                }) 
            }
    })
    

    
    output$seqPlot <- renderPlot({ 
        print("about to make plot")
        obj <- orderObject
        # print(str(as.list(isolate(coordinatesObject))))
        makePlot(obj,isolate(coordinatesObject))
        print("done making plot")
        })
    
    output$down <- downloadHandler(
        filename = function(){
            "plot.png"
        },
        content = function(file){
            if (input$filetype == "PNG") png(file)
            if (input$filetype == "SVG") svg(file)
            print("download handler")
            makePlot(orderObject, coordinatesObject, plotFAST = FALSE)
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
#' @import shiny
#' @export
methylScaper <- function() {shinyApp(ui = ui, server = server)}
