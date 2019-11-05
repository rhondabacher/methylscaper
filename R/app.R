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
           downloadButton("down", label = "Download the plot")
          
        )
    )
)

server <- function(input, output) {
    output$seqPlot <- renderPlot({
        
        input.GCH <- day7$gch
        input.HCG <- day7$hcg
        
        if (!is.null(input$gch.file) & !is.null(input$hcg.file))
        {
            input.GCH <- read.table(input$gch.file$datapath, header=T, row.names = 1, stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
            input.HCG <- read.table(input$hcg.file$datapath, header=T, row.names = 1, stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
            
        }
        
        processed.brush <- list(first.row = 0, last.row = 0, first.col = 0, last.col = 0, weight.color = "red")
        if (!is.null(input$plot_brush))
        {
            n <- nrow(input.GCH)
            m <- ncol(input.GCH)
            processed.brush <- handleBrushCoordinates(input$plot_brush,input$brush.choice, n, m)
        }
        orderObj <- buildOrderObjectShiny(input.GCH, input.HCG, method = input$method,
                                          weight.start = processed.brush$first.col, 
                                          weight.stop = processed.brush$last.col, 
                                          weight.color = processed.brush$weight.color)
        refineOrderShiny(orderObj,
                         refine.method = input$refineMethod)
        makePlot(orderObj)
        })
    
    output$down <- downloadHandler(
        filename = function(){
            "test.png"
        },
        content = function(file){
            if (input$filetype == "PNG") png(file)
            if (input$filetype == "SVG") svg(file)
            
            input.GCH <- day7$gch
            input.HCG <- day7$hcg
            
            if (!is.null(input$gch.file) & !is.null(input$hcg.file))
            {
                input.GCH <- read.table(input$gch.file$datapath, header=T, row.names = 1, stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
                input.HCG <- read.table(input$hcg.file$datapath, header=T, row.names = 1, stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
                
            } 
        
            processed.brush <- list(first.row = 0, last.row = 0, first.col = 0, last.col = 0, weight.color = "red")
            if (!is.null(input$plot_brush))
            {
                n <- nrow(input.GCH)
                m <- ncol(input.GCH)
                processed.brush <- handleBrushCoordinates(input$plot_brush, input$brush.choice, n, m)
            }
            
            orderObj <- buildOrderObjectShiny(input.GCH, input.HCG, method = input$method,
                                              weight.start = processed.brush$first.col, 
                                              weight.stop = processed.brush$last.col, 
                                              weight.color = processed.brush$weight.color)
            refineOrderShiny(orderObj,
                             refine.method = input$refineMethod)
            makePlot(orderObj, plotFAST = FALSE)
            dev.off()
        }
    )
    
    output$info <- renderText({
        input.GCH <- day7$gch
        
        if (!is.null(input$gch.file) & !is.null(input$hcg.file))
        {
            input.GCH <- read.table(input$gch.file$datapath, header=T, row.names = 1, stringsAsFactors = F, quote = "", sep = "\t", comment.char = "")
            
        }
        processed.brush <- list(first.row = 0, last.row = 0, first.col = 0, last.col = 0)
        
        if (!is.null(input$plot_brush))
        {
            n <- nrow(input.GCH)
            m <- ncol(input.GCH)
            processed.brush <- handleBrushCoordinates(input$plot_brush, input$brush.choice, n, m)
        }
        
        paste0("Refinement selection: ", refine.start.g, " ", refine.stop.g, "\n",
               "Weighting selection: ", weight.start.g, " ", weight.stop.g)
    })
}





# Run the application 
#' @import shiny
#' @export
methylScaper <- function() {shinyApp(ui = ui, server = server)}
