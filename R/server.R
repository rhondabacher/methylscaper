server <- function(input, output) {

  actionsLog <- reactiveValues(log = c("")) # logs the actions taken wrt the plot


  ####################
  # Single-cell data #
  ####################

  sc_seq_data <- reactiveValues(gch = NULL, hcg = NULL) # for raw data
  sc_input_data <- reactiveValues(gch = NULL, hcg = NULL) # for state matrices

  sc_input_folder <- reactiveValues(path = NULL)
  

  ## preprocessing tab
  observe({
   volumes = getVolumes()
    shinyDirChoose(input, 'folder', roots=volumes())
    path.list <- input$folder["path"][[1]]
    path <- paste(unlist(path.list), collapse = "/")
    print(path)
    if (!is.null(path.list)) sc_input_folder$path <- path

   })
  
  output$sc_folder_name <- renderText({
    if (is.null(sc_input_folder$path)) "Choose an input folder."
    else sc_input_folder$path
    })

  observeEvent(input$run.subset,{
   if (!is.na(sc_input_folder$path))
     {
       progress <- Progress$new()
       progress$set(message = "Loading single cell data", value = 0)
       on.exit(progress$close())

       updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
            progress$set(value = value, message = message, detail = detail)}

       print(paste("Begin SC processing in", sc_input_folder$path, "at chr", input$chromosome.number))
       dat.subset <- subsetSC(sc_input_folder$path, input$chromosome.number, updateProgress = updateProgress) # this is very slow
       print("Done with subset, beginning prepSC")
       prepsc.out <- prepSC(dat.subset$gc.seq.sub, dat.subset$cg.seq.sub,
                            startPos = 105636488, endPos = 105636993, updateProgress = updateProgress)
       print("Done with single cell processing")
     }
  })

  ## seriation tab
    
  observe({
    if (!is.null(input$gch_seq_file) & !is.null(input$hcg_seq_file))
    {
      isolate({
        progress <- Progress$new()
        progress$set(message = "Loading GC data", value = 0)
        on.exit(progress$close())

        sc_seq_data$gch <- readRDS(input$gch_seq_file$datapath)
        progress$set(message = "Loading CG data", value = 0.5)
        sc_seq_data$hcg <- readRDS(input$hcg_seq_file$datapath)
        actionsLog$log <- c(actionsLog$log, paste("Loading GCH RDS file:", input$gch_seq_file$name))
        actionsLog$log <- c(actionsLog$log, paste("Loading HCG RDS file:", input$hcg_seq_file$name))
      })
    }
  })

  output$startPos <- renderUI({
      if (!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg))
      {
          cg.max.pos <- max(sapply(sc_seq_data$hcg, FUN=function(x) max(x$pos)))
          cg.min.pos <- min(sapply(sc_seq_data$hcg, FUN=function(x) min(x$pos)))
          gc.max.pos <- max(sapply(sc_seq_data$gch, FUN=function(x) max(x$pos)))
          gc.min.pos <- min(sapply(sc_seq_data$gch, FUN=function(x) min(x$pos)))

          start <- pmax(cg.min.pos, gc.min.pos)
          numericInput(inputId = "startPos", label = "Start Position", min = 0, 
                      value = start)
      }

  })
  
  output$endPos <- renderUI({
      if (!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg))
      {
          cg.max.pos <- max(sapply(sc_seq_data$hcg, FUN=function(x) max(x$pos)))
          cg.min.pos <- min(sapply(sc_seq_data$hcg, FUN=function(x) min(x$pos)))
          gc.max.pos <- max(sapply(sc_seq_data$gch, FUN=function(x) max(x$pos)))
          gc.min.pos <- min(sapply(sc_seq_data$gch, FUN=function(x) min(x$pos)))

          end <- pmax(cg.min.pos, gc.min.pos) + 1000
          numericInput(inputId = "endPos", label = "End Position", min = 0, 
                      value = end)
      }

  })
    output$positionSlider <- renderUI({
        if (!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg))
        {
            cg.max.pos <- max(sapply(sc_seq_data$hcg, FUN=function(x) max(x$pos)))
            cg.min.pos <- min(sapply(sc_seq_data$hcg, FUN=function(x) min(x$pos)))
            gc.max.pos <- max(sapply(sc_seq_data$gch, FUN=function(x) max(x$pos)))
            gc.min.pos <- min(sapply(sc_seq_data$gch, FUN=function(x) min(x$pos)))
            start <- input$startPos
            end <- input$endPos
            
            if (start < cg.min.pos | start < gc.min.pos | end > cg.max.pos | end > gc.max.pos)
            {
                showNotification("Selected range is out of bounds. Please choose a valid starting and end position to generate the plot.", type="error")
                return(NULL)
            }
            if (end -  start > 10000)
            {
                showNotification("Selected range is longer than 10k bp, reducing length for stability.", type="warning")
                end <- start + 10000
            }
            len <- end - start

            sliderInput(inputId = "positionSliderInput", label = "Position adjustment slider", min = start - len, max = end + len,
                        value = c(start, end))
        }

    })

  observe({
    if (!is.null(input$positionSliderInput))
    {
        progress <- Progress$new()
        progress$set(message = "Beginning single-cell processing", value = 0)
        on.exit(progress$close())

        updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
            progress$set(value = value, message = message, detail = detail)}

        prep_out <- prepSC(sc_seq_data$gch, sc_seq_data$hcg, input$positionSliderInput[1],
                           input$positionSliderInput[2],
                           updateProgress = updateProgress)
        if (!is.list(prep_out)) {
          showNotification("No valid sites in designated range. Try different start and end positions or a larger range.")
         } else {
         temp.gch <- prep_out$gch
         temp.hcg <- prep_out$hcg
          if (nrow(temp.gch) == nrow(temp.hcg))
          {
            sc_coordinatesObject$refine.start <- 0
            sc_coordinatesObject$refine.stop <- 0
            sc_coordinatesObject$weight.start <- 0
            sc_coordinatesObject$weight.stop <- 0
            sc_input_data$gch <- temp.gch
            sc_input_data$hcg <- temp.hcg
            isolate({
              actionsLog$log <- c(actionsLog$log, paste("Beginning single-cell data analysis"))
              actionsLog$log <- c(actionsLog$log, paste("From position", input$positionSliderInput[1], "to", input$positionSliderInput[2]))
            })
          }
         }


    }

  })

  # this object keeps track of the coordinates for refinement and weighting
  sc_coordinatesObject <- reactiveValues(refine.start = 0, refine.stop = 0,
                                         weight.start = 0, weight.stop = 0, weight.color = "red")
  # now construct the sc_orderObject
  sc_orderObject <- reactiveValues(toClust = 0, order1 = 0)
  observe({ if (!is.null(sc_input_data$gch) & !is.null(sc_input_data$hcg))
  {
    progress <- Progress$new()
    progress$set(message = "Beginning seriation", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      progress$set(value = value, message = message, detail = detail)}

    tempObj <- buildOrderObjectShiny(sc_input_data$gch, sc_input_data$hcg, input$sc_ser_method, sc_coordinatesObject, updateProgress)
    sc_orderObject$order1 <- tempObj$order1
    sc_orderObject$toClust <- tempObj$toClust
    isolate({
      actionsLog$log <- c(actionsLog$log,
                          paste("Ordering with", input$sc_ser_method))
    })
  }

  })

  # this handles updates to sc_coordinatesObject
  observeEvent(input$sc_plot_brush, {
    n <- nrow(sc_input_data$gch)
    m <- ncol(sc_input_data$hcg)
    processed.brush <- handleBrushCoordinates(input$sc_plot_brush, n, m)

    if (isolate(input$sc_brush_choice) == "Weighting")
    {
      sc_coordinatesObject$refine.start <- 0
      sc_coordinatesObject$refine.stop <- 0
      sc_coordinatesObject$weight.start <- processed.brush$first.col
      sc_coordinatesObject$weight.stop <- processed.brush$last.col
      sc_coordinatesObject$weight.color <- processed.brush$weight.color
      isolate({
        actionsLog$log <- c(actionsLog$log,
                            paste("Weighting", processed.brush$weight.color, "columns",
                                  processed.brush$first.col, "to",
                                  processed.brush$last.col))
      })
    }
    if (isolate(input$sc_brush_choice) == "Refinement")
    {
      s <- processed.brush$first.row
      f <- processed.brush$last.row
      if (s != f)
      {
        sc_coordinatesObject$refine.start <- s
        sc_coordinatesObject$refine.stop <- f
        sc_orderObject$order1 <- refineOrderShiny(isolate(sc_orderObject),
                                                  refine.method = isolate(input$sc_refine_method),
                                                  sc_coordinatesObject)
        isolate({
          actionsLog$log <- c(actionsLog$log,
                              paste("Refining rows",
                                    processed.brush$first.row, "to",
                                    processed.brush$last.row))
          actionsLog$log <- c(actionsLog$log,
                              paste("Applying refinement with", input$sc_refine_method))
        })
      }
    }


  })

  observeEvent( input$sc_force_reverse, {
    isolate({
      if (sc_coordinatesObject$refine.start == sc_coordinatesObject$refine.stop)
      {
        sc_orderObject$order1 <- rev(sc_orderObject$order1)
        actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to", nrow(sc_input_data$gch)))
      }
      else
      {
        sc_orderObject$order1[sc_coordinatesObject$refine.start : sc_coordinatesObject$refine.stop] <-
          sc_orderObject$order1[sc_coordinatesObject$refine.stop : sc_coordinatesObject$refine.start]
        actionsLog$log <- c(actionsLog$log,
                            paste("Reversing rows", sc_coordinatesObject$refine.start,
                                  "to", sc_coordinatesObject$refine.stop))

      }
    })
  })



  output$sc_seqPlot <- renderPlot({
    obj <- sc_orderObject
    if (sum(obj$toClust) == 0) {showNotification("Select methylation data files to generate the plot.", type="message");NULL}
    else makePlot(obj,isolate(sc_coordinatesObject))
  }, height=600, width=600)

  output$sc_plot_down <- downloadHandler(
    filename = function(){
      if (input$sc_filetype == "PNG") return("plot.png")
      if (input$sc_filetype == "SVG") return("plot.svg")
      if (input$sc_filetype == "PDF") return("plot.pdf")
    },
    content = function(file){
      if (input$sc_filetype == "PNG") png(file)
      if (input$sc_filetype == "SVG") svglite::svglite(file)
      if (input$sc_filetype == "PDF") pdf(file)

      makePlot(sc_orderObject, sc_coordinatesObject, drawLines = FALSE, plotFAST = FALSE)
      dev.off()
    }
  )

  output$sc_log_down <- downloadHandler(
    filename = function(){
      "changes.txt"
    },
    content = function(file){
      fileConn <- file(file)
      writeLines(actionsLog$log, fileConn)
      close(fileConn)
    }
  )

  output$sc_info <- renderText({
    paste0("Refinement selection: ", sc_coordinatesObject$refine.start, " ", sc_coordinatesObject$refine.stop, "\n",
           "Weighting selection: ", sc_coordinatesObject$weight.start, " ", sc_coordinatesObject$weight.stop)
  })

  output$sc_proportion_color_histogram <- renderPlot({
    obj <- sc_orderObject
    if (sum(obj$toClust) == 0)
    {showNotification("Select methylation data files to generate the plot.", type="message");NULL}
    else proportion_color(obj, plotHistogram = TRUE, color = toupper(input$sc_proportion_choice))
  })

  output$sc_proportion_hist_download <- downloadHandler(
    filename = function(){
      if (input$filetype == "PNG") return("hist.png")
      if (input$filetype == "SVG") return("hist.svg")
      if (input$filetype == "PDF") return("hist.pdf")
    },
    content = function(file){
      if (input$filetype == "PNG") png(file)
      if (input$filetype == "SVG") svglite::svglite(file)
      if (input$filetype == "PDF") pdf(file)

      proportion_color(sc_orderObject, plotHistogram = TRUE,
                       color = toupper(input$sc_proportion_choice))
      dev.off()
    }
  )
  output$sc_proportion_data_download <- downloadHandler(
    filename = function(){
      return("proportion_data.csv")
    },
    content = function(file){
      dat <-  proportion_color(sc_orderObject, plotHistogram = FALSE,
                               color = toupper(input$sc_proportion_choice))
      write.csv(dat, file = file)
    }
  )

  output$sc_percent_C <- renderPlot({
    obj <- sc_orderObject
    if (sum(obj$toClust) == 0)
    {showNotification("Select methylation data files to generate the plot.", type="message");NULL}
    else percent_C(obj, plotPercents=TRUE)
  })

  output$sc_percentC_plot_download <- downloadHandler(
    filename = function(){
      if (input$filetype == "PNG") return("percentC.png")
      if (input$filetype == "SVG") return("percentC.svg")
      if (input$filetype == "PDF") return("percentC.pdf")
    },
    content = function(file){
      if (input$filetype == "PNG") png(file)
      if (input$filetype == "SVG") svglite::svglite(file)
      if (input$filetype == "PDF") pdf(file)

      percent_C(sc_orderObject, plotPercents = TRUE)
      dev.off()
    }
  )

  output$sc_percentC_data_download <- downloadHandler(
    filename = function(){
      return("proportion_data.RData")
    },
    content = function(file){
      dat <-  percent_C(sc_orderObject, plotPercents = FALSE)
      save(dat, file = file)
    }
  )




  ########################
  # Single-molecule data #
  ########################
  sm_input_data <- reactiveValues(gch = NULL, hcg = NULL)

  # alignment handling
  observeEvent(input$run.align, {
    ref <- read.fasta(input$ref.file$datapath)
    fasta <- read.fasta(input$fasta.file$datapath)

    progress <- Progress$new()
    progress$set(message = "Beginning alignment", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      progress$set(value = value, message = message, detail = detail)}

    align.out <- runAlign(ref, fasta, updateProgress = updateProgress,
                          log.file = input$processing.log.name)

    hcg.file.name <- input$hcg.file.name
    gch.file.name <- input$gch.file.name

    writeMethylationData(dat = align.out$hcg, filepath = hcg.file.name)
    writeMethylationData(dat = align.out$gch, filepath = gch.file.name)

  })

  observe({if (!is.null(input$sm_gch_file) & !is.null(input$sm_hcg_file))
  {
    temp.gch <- readMethylationData(filepath = input$sm_gch_file$datapath)
    temp.hcg <- readMethylationData(filepath = input$sm_hcg_file$datapath)
    if (nrow(temp.gch) == nrow(temp.hcg))
    {
      sm_coordinatesObject$refine.start <- 0
      sm_coordinatesObject$refine.stop <- 0
      sm_coordinatesObject$weight.start <- 0
      sm_coordinatesObject$weight.stop <- 0
      sm_input_data$gch <- temp.gch
      sm_input_data$hcg <- temp.hcg
      sm_input_data$datatype <- "sm"
      isolate({
        actionsLog$log <- c(actionsLog$log, paste("Beginning single-molecule data analysis"))
        actionsLog$log <- c(actionsLog$log, paste("Loading GCH file:", input$gch.file$name))
        actionsLog$log <- c(actionsLog$log, paste("Loading HCG file:", input$hcg.file$name))
      })
    }

  }})

  # this object keeps track of the coordinates for refinement and weighting
  sm_coordinatesObject <- reactiveValues(refine.start = 0, refine.stop = 0,
                                      weight.start = 0, weight.stop = 0, weight.color = "red")
  # now construct the sm_orderObject
  sm_orderObject <- reactiveValues(toClust = 0, order1 = 0)
  observe({ if (!is.null(sm_input_data$gch) & !is.null(sm_input_data$hcg))
  {
    progress <- Progress$new()
    progress$set(message = "Beginning seriation", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      progress$set(value = value, message = message, detail = detail)}

    tempObj <- buildOrderObjectShiny(sm_input_data$gch, sm_input_data$hcg, input$sm_ser_method, sm_coordinatesObject, updateProgress)
    sm_orderObject$order1 <- tempObj$order1
    sm_orderObject$toClust <- tempObj$toClust
    isolate({
      actionsLog$log <- c(actionsLog$log,
                          paste("Ordering with", input$sm_ser_method))
    })
  }

  })

  # this handles updates to sm_coordinatesObject
  observeEvent(input$sm_plot_brush, {
    n <- nrow(sm_input_data$gch)
    m <- ncol(sm_input_data$hcg)
    processed.brush <- handleBrushCoordinates(input$sm_plot_brush, n, m)

    if (isolate(input$sm_brush_choice) == "Weighting")
    {
      sm_coordinatesObject$refine.start <- 0
      sm_coordinatesObject$refine.stop <- 0
      sm_coordinatesObject$weight.start <- processed.brush$first.col
      sm_coordinatesObject$weight.stop <- processed.brush$last.col
      sm_coordinatesObject$weight.color <- processed.brush$weight.color
      isolate({
        actionsLog$log <- c(actionsLog$log,
                            paste("Weighting", processed.brush$weight.color, "columns",
                                  processed.brush$first.col, "to",
                                  processed.brush$last.col))
      })
    }
    if (isolate(input$sm_brush_choice) == "Refinement")
    {
      s <- processed.brush$first.row
      f <- processed.brush$last.row
      if (s != f)
      {
        sm_coordinatesObject$refine.start <- s
        sm_coordinatesObject$refine.stop <- f
        sm_orderObject$order1 <- refineOrderShiny(isolate(sm_orderObject),
                                               refine.method = isolate(input$sm_refine_method),
                                               sm_coordinatesObject)
        isolate({
          actionsLog$log <- c(actionsLog$log,
                              paste("Refining rows",
                                    processed.brush$first.row, "to",
                                    processed.brush$last.row))
          actionsLog$log <- c(actionsLog$log,
                              paste("Applying refinement with", input$sm_refine_method))
        })
      }
    }


  })

  observeEvent( input$sm_force_reverse, {
    isolate({
      if (sm_coordinatesObject$refine.start == sm_coordinatesObject$refine.stop)
      {
        sm_orderObject$order1 <- rev(sm_orderObject$order1)
        actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to", nrow(sm_input_data$gch)))
      }
      else
      {
        sm_orderObject$order1[sm_coordinatesObject$refine.start : sm_coordinatesObject$refine.stop] <-
          sm_orderObject$order1[sm_coordinatesObject$refine.stop : sm_coordinatesObject$refine.start]
        actionsLog$log <- c(actionsLog$log,
                            paste("Reversing rows", sm_coordinatesObject$refine.start,
                                  "to", sm_coordinatesObject$refine.stop))

      }
    })
  })



  output$sm_seqPlot <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0) {showNotification("Select methylation data files to generate the plot.", type="message");NULL}
    else makePlot(obj,isolate(sm_coordinatesObject))
  }, height=600, width=600)

  output$sm_plot_down <- downloadHandler(
    filename = function(){
      if (input$sm_filetype == "PNG") return("plot.png")
      if (input$sm_filetype == "SVG") return("plot.svg")
      if (input$sm_filetype == "PDF") return("plot.pdf")
    },
    content = function(file){
      if (input$sm_filetype == "PNG") png(file)
      if (input$sm_filetype == "SVG") svglite::svglite(file)
      if (input$sm_filetype == "PDF") pdf(file)

      makePlot(sm_orderObject, sm_coordinatesObject, drawLines = FALSE, plotFAST = FALSE)
      dev.off()
    }
  )

  output$sm_log_down <- downloadHandler(
    filename = function(){
      "changes.txt"
    },
    content = function(file){
      fileConn <- file(file)
      writeLines(actionsLog$log, fileConn)
      close(fileConn)
    }
  )

  output$sm_info <- renderText({
    paste0("Refinement selection: ", sm_coordinatesObject$refine.start, " ", sm_coordinatesObject$refine.stop, "\n",
           "Weighting selection: ", sm_coordinatesObject$weight.start, " ", sm_coordinatesObject$weight.stop)
  })

  output$sm_proportion_color_histogram <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0)
    {showNotification("Select methylation data files to generate the plot.", type="message");NULL}
    else proportion_color(obj, plotHistogram = TRUE, color = toupper(input$sm_proportion_choice))
  })

  output$sm_proportion_hist_download <- downloadHandler(
    filename = function(){
      if (input$filetype == "PNG") return("hist.png")
      if (input$filetype == "SVG") return("hist.svg")
      if (input$filetype == "PDF") return("hist.pdf")
    },
    content = function(file){
      if (input$filetype == "PNG") png(file)
      if (input$filetype == "SVG") svglite::svglite(file)
      if (input$filetype == "PDF") pdf(file)

      proportion_color(sm_orderObject, plotHistogram = TRUE,
                       color = toupper(input$sm_proportion_choice))
      dev.off()
    }
  )
  output$sm_proportion_data_download <- downloadHandler(
    filename = function(){
      return("proportion_data.csv")
    },
    content = function(file){
      dat <-  proportion_color(sm_orderObject, plotHistogram = FALSE,
                               color = toupper(input$sm_proportion_choice))
      write.csv(dat, file = file)
    }
  )

  output$sm_percent_C <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0)
    {showNotification("Select methylation data files to generate the plot.", type="message");NULL}
    else percent_C(obj, plotPercents=TRUE)
  })

  output$sm_percentC_plot_download <- downloadHandler(
    filename = function(){
      if (input$filetype == "PNG") return("percentC.png")
      if (input$filetype == "SVG") return("percentC.svg")
      if (input$filetype == "PDF") return("percentC.pdf")
    },
    content = function(file){
      if (input$filetype == "PNG") png(file)
      if (input$filetype == "SVG") svglite::svglite(file)
      if (input$filetype == "PDF") pdf(file)

      percent_C(sm_orderObject, plotPercents = TRUE)
      dev.off()
    }
  )

  output$sm_percentC_data_download <- downloadHandler(
    filename = function(){
      return("proportion_data.RData")
    },
    content = function(file){
      dat <-  percent_C(sm_orderObject, plotPercents = FALSE)
      save(dat, file = file)
    }
  )
}
