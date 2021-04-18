server <- function(input, output, session) {
		
  actionsLog <- reactiveValues(log = c("")) # logs the actions taken wrt the plot



  ####################
  # Single-cell data #
  ####################

  sc_seq_data <- reactiveValues(gch = NULL, hcg = NULL) # for raw data
  sc_raw_data <- reactiveValues(gch = NULL, hcg = NULL)
  sc_input_data <- reactiveValues(gch = NULL, hcg = NULL) # for state matrices
  sc_input_folder <- reactiveValues(path = NULL)
  mouse_bm <- NULL
  human_bm <- NULL
  singlecell_subset <- NULL
	singlemolecule_example <- NULL
	outname <- reactiveValues(usename = "example_data")
	outname_rds <- reactiveValues(usename = "example_data")
	
  ## preprocessing tab
  observe({
    if (is.null(input$sc_rds_file) & input$seriate_sc == "Preprocessing"  & input$big_tab == "Single-cell")
    {
      showNotification("Please select the input files to begin", 
                                      type="message", duration=4)
    }
   })

  observeEvent(input$run_subset,{
		
    validate(need(!is.null(input$sc_met_files$name[1]) & !is.null(input$sc_acc_files$name[1]),
    							message = "Please choose an input directory.", label = "sc_input_folder"))

    progress <- Progress$new()
    progress$set(message = "Loading single-cell data", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
         progress$set(value = value, message = message, detail = detail) }
    showNotification(paste("Begin SC processing at chr", input$chromosome_number))
		
    dat_subset <- subsetSC(list(input$sc_met_files$datapath, input$sc_acc_files$datapath,
																input$sc_met_files$data, input$sc_met_files$data), 
													input$chromosome_number, updateProgress = updateProgress) 

    showNotification("Done with single cell processing")
    sc_raw_data$gch <- dat_subset$gch
    sc_raw_data$hcg <- dat_subset$hcg
    rm(dat_subset)
    showNotification("Removed temporary raw data; Click button to download now.", duration=3)
  })
	
	
	observe({
	    if (is.null(sc_raw_data$gch)){
				shinyjs::disable("sc_preprocessing_down") 
	   } else {
	      shinyjs::enable("sc_preprocessing_down")
			}
	  })
		
output$sc_preprocessing_down <- downloadHandler(
  filename = function(){
      paste0("methylscaper_preprocessed_",input$chromosome_number,".rds")
    },
    content = function(file){
      validate(need(!is.null(sc_raw_data$gch) & !is.null(sc_raw_data$hcg), 
									message = "Data has not been processed."))
      print("Saving data")
      saveRDS(list(gch = sc_raw_data$gch, hcg = sc_raw_data$hcg), file = file)
      sc_raw_data$gch <- NULL
      sc_raw_data$hcg <- NULL
    }
  )
 
  ## Visualization tab
    
 
  observe({
      if (is.null(input$sc_rds_file) & input$seriate_sc == "Visualization" & input$big_tab == "Single-cell")

      {
        showNotification("Provide select the RDS file to begin", 
                                        type="message", duration=4)
      }
    if (!is.null(input$sc_rds_file))
    {
      isolate({
        progress <- Progress$new()
        progress$set(message = "Loading data", value = 0)
        on.exit(progress$close())
        temp <- readRDS(input$sc_rds_file$datapath)
        sc_seq_data$gch <- temp$gch
        sc_seq_data$hcg <- temp$hcg
				outname_rds$usename <- tools::file_path_sans_ext(input$sc_rds_file$name)
        actionsLog$log <- c(actionsLog$log, paste("Loading data:",
                                                input$sc_rds_file$name))
      })
			if(is.null(input$organism_choice)) {
      showNotification("Now select Organism and begin selecting genes", 
                                      type="message", duration=10)
																		}
    }
  })

  # Genes <- reactiveValues()
  observeEvent(input$organism_choice, {
      
      if(!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg)) {
      
        if (input$organism_choice == "Human") {
            data("human_bm", package="methylscaper", envir = environment())
            getchr <- sc_seq_data$gch[[1]]$chr[1]
            cg_max_pos <- suppressWarnings(max(vapply(sc_seq_data$hcg, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            cg_min_pos <- suppressWarnings(min(vapply(sc_seq_data$hcg, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_max_pos <- suppressWarnings(max(vapply(sc_seq_data$gch, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_min_pos <- suppressWarnings(min(vapply(sc_seq_data$gch, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))
			getmin <- pmin(cg_min_pos, gc_min_pos) 
			getmin <- max(c(0,getmin - 100000))
			getmax <- pmax(cg_max_pos, gc_max_pos) 
			getmax <- max(c(0,getmax + 100000))
            hum_bm_sub <- subset(human_bm, 
								human_bm$chromosome_name == getchr & 
								human_bm$start_position >= getmin & 
								human_bm$end_position <= getmax)
             Genes <- sort(unique(hum_bm_sub$hgnc_symbol))
        } else if (input$organism_choice == "Mouse") {
            data("mouse_bm", package="methylscaper", envir = environment())
            getchr <- sc_seq_data$gch[[1]]$chr[1]
            cg_max_pos <- suppressWarnings(max(vapply(sc_seq_data$hcg, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            cg_min_pos <- suppressWarnings(min(vapply(sc_seq_data$hcg, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_max_pos <- suppressWarnings(max(vapply(sc_seq_data$gch, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_min_pos <- suppressWarnings(min(vapply(sc_seq_data$gch, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))
			getmin <- pmin(cg_min_pos, gc_min_pos) 
			getmin <- max(c(0,getmin - 100000))
			getmax <- pmax(cg_max_pos, gc_max_pos) 
			getmax <- max(c(0,getmax + 100000))

            mouse_bm_sub <- subset(mouse_bm, mouse_bm$chromosome_name == getchr & 
								mouse_bm$start_position >= getmin & 
								mouse_bm$end_position <= getmax)
            Genes <- sort(unique(mouse_bm_sub$mgi_symbol))
        } else if (input$organism_choice == "Other") {
            Genes = "Click here to begin manual start and end selection."
        }
     updateSelectizeInput(session, "geneList",
                               choices = Genes,
                               server = TRUE, selected = ' ')
     }
   })
   
   
   
  output$startPos <- renderUI({
     
      if (!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg) & input$geneList != "") {
          
          if (input$organism_choice == "Mouse") {
	          data("mouse_bm", package="methylscaper", envir = environment())
              gene_select <- subset(mouse_bm, mouse_bm$mgi_symbol == input$geneList)
          }
          if (input$organism_choice == "Human") {
	          data("human_bm", package="methylscaper", envir = environment())
              gene_select <- subset(human_bm, human_bm$hgnc_symbol == input$geneList)
          }
          if (input$organism_choice == "Other") {
            cg_max_pos <- suppressWarnings(max(vapply(sc_seq_data$hcg, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            cg_min_pos <- suppressWarnings(min(vapply(sc_seq_data$hcg, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_max_pos <- suppressWarnings(max(vapply(sc_seq_data$gch, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_min_pos <- suppressWarnings(min(vapply(sc_seq_data$gch, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))


            start <- pmax(cg_min_pos, gc_min_pos)
            gene_select <- data.frame(start_position = start)
          }
          start <- gene_select$start_position
          numericInput(inputId = "startPos", label = "Start Position", min = 0,
                      value = start)
      }

  })
  
  output$endPos <- renderUI({
    
        if (!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg) & input$geneList != "") {

          if (input$organism_choice == "Mouse") {
              data("mouse_bm", package="methylscaper", envir = environment())
              gene_select <- subset(mouse_bm, mouse_bm$mgi_symbol == input$geneList)
          }
          if (input$organism_choice == "Human") {
              data("human_bm", package="methylscaper", envir = environment())
	      	  gene_select <- subset(human_bm, human_bm$hgnc_symbol == input$geneList)
          }
          if (input$organism_choice == "Other") {
            cg_max_pos <- suppressWarnings(max(vapply(sc_seq_data$hcg, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            cg_min_pos <- suppressWarnings(min(vapply(sc_seq_data$hcg, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_max_pos <- suppressWarnings(max(vapply(sc_seq_data$gch, FUN=function(x) {max(x$pos, na.rm=TRUE)}, numeric(1))))
            gc_min_pos <- suppressWarnings(min(vapply(sc_seq_data$gch, FUN=function(x) {min(x$pos, na.rm=TRUE)}, numeric(1))))

            end <- pmax(cg_min_pos, gc_min_pos) + 5000
            gene_select <- data.frame(end_position = end)
         }
          end <- gene_select$end_position
          numericInput(inputId = "endPos", label = "End Position", min = 0, 
                      value = end)
      }

  })
    output$positionSlider <- renderUI({
        if (!is.null(sc_seq_data$gch) & !is.null(sc_seq_data$hcg) & input$geneList != "") {
		    isolate({
		      actionsLog$log <- c(actionsLog$log,
		                          paste("Current gene selected: ", input$geneList))
		    })

		if (!is.null(input$startPos) & !is.null(input$endPos)) {
			start <- input$startPos
			end <- input$endPos
            
        if (end -  start > 50000) {
            showNotification("Selected range is longer than 50k bp, plot may take a few 
                        seconds to render", duration=3)
        }
        if (end -  start > 100000) {
            showNotification("Selected range is longer than 100k bp, this is not optimal for 
                        visualization, reducing to 100k bp.", duration=10)

            end <- start + 100000
        }
        if (start > end) {
            end <- start + 2000
        }
        len <- end - start
		if (len > 0) {
        sliderInput(inputId = "positionSliderInput", 
                    label = "Position adjustment slider", 
                    min = start - len, max = end + len,
                        value = c(start, end))
					}
		}
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


		
				prep_out <- prepSC(sc_seq_data, 
                            input$positionSliderInput[1],
                            input$positionSliderInput[2],
                            updateProgress = updateProgress)
        if (!is.list(prep_out)) {

					showElement("sc_seqPlot")
					toggleElement("sc_seqPlot")
          showNotification("No valid sites in designated range. Choose another gene or adjust  
                      start and end positions with a larger range.", duration=5)
		    isolate({
		      actionsLog$log <- c(actionsLog$log,
		                          paste("No valid sites for gene ", input$geneList))
		    })
         
			 } else {
				 showElement("sc_seqPlot")
				 temp_gch <- prep_out$gch

         temp_hcg <- prep_out$hcg
         if (nrow(temp_gch) == nrow(temp_hcg)) {
            sc_coordinatesObject$refine_start <- 0
            sc_coordinatesObject$refine_stop <- 0
            sc_coordinatesObject$weight_start <- 0
            sc_coordinatesObject$weight_stop <- 0
            sc_input_data$gch <- temp_gch
            sc_input_data$hcg <- temp_hcg
            isolate({
              actionsLog$log <- c(actionsLog$log, paste("Beginning single-cell data analysis"))
              actionsLog$log <- c(actionsLog$log, paste("From position",
                               input$positionSliderInput[1], 
                               "to", input$positionSliderInput[2]))
            })
          }
         }
    }

  })
	

# this object keeps track of the coordinates for refinement and weighting
  sc_coordinatesObject <- reactiveValues(refine_start = 0, refine_stop = 0,
                                         weight_start = 0, weight_stop = 0, 
                                         weight_color = "red")
  # now construct the sc_orderObject
  sc_orderObject <- reactiveValues(toClust = 0, order1 = 0)
  observe({ if (!is.null(sc_input_data$gch) & !is.null(sc_input_data$hcg))
  {
    progress <- Progress$new()
    progress$set(message = "Beginning seriation", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      progress$set(value = value, message = message, detail = detail)}

    tempObj <- buildOrderObjectShiny(sc_input_data,
                         input$sc_ser_method, sc_coordinatesObject, 
                         updateProgress)
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
    validate(need(sc_input_data$gch, "Please provide input data"))
    n <- nrow(sc_input_data$gch)
    validate(need(sc_input_data$hcg, "Please provide input data"))  
    m <- ncol(sc_input_data$hcg)
    processed_brush <- handleBrushCoordinates(input$sc_plot_brush, n, m)

    if (isolate(input$sc_brush_choice) == "Weighting")
    {
      sc_coordinatesObject$refine_start <- 0
      sc_coordinatesObject$refine_stop <- 0
      sc_coordinatesObject$weight_start <- processed_brush$first_col
      sc_coordinatesObject$weight_stop <- processed_brush$last_col
      sc_coordinatesObject$weight_color <- processed_brush$weight_color
      isolate({
        actionsLog$log <- c(actionsLog$log,
                            paste("Weighting", processed_brush$weight_color,
                                  "columns",
                                  processed_brush$first_col, "to",
                                  processed_brush$last_col))
      })
    }
    if (isolate(input$sc_brush_choice) == "Refinement")
    {
      s <- processed_brush$first_row
      f <- processed_brush$last_row
      if (s != f)
      {
        sc_coordinatesObject$refine_start <- s
        sc_coordinatesObject$refine_stop <- f
        sc_orderObject$order1 <- refineOrderShiny(isolate(sc_orderObject),
                                        refine_method = isolate(input$sc_refine_method),
                                        sc_coordinatesObject)
        isolate({
          actionsLog$log <- c(actionsLog$log,
                              paste("Refining rows",
                                    processed_brush$first_row, "to",
                                    processed_brush$last_row))
          actionsLog$log <- c(actionsLog$log,
                              paste("Applying refinement with", 
                              input$sc_refine_method))
        })
      }
    }


  })

  observeEvent( input$sc_force_reverse, {
    isolate({
      if (sc_coordinatesObject$refine_start == sc_coordinatesObject$refine_stop)
      {
        sc_orderObject$order1 <- rev(sc_orderObject$order1)
        actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to", 
                                                nrow(sc_input_data$gch)))
      }
      else
      {
        sc_orderObject$order1[sc_coordinatesObject$refine_start : sc_coordinatesObject$refine_stop] <-
          sc_orderObject$order1[sc_coordinatesObject$refine_stop : sc_coordinatesObject$refine_start]
        actionsLog$log <- c(actionsLog$log,
                            paste("Reversing rows", sc_coordinatesObject$refine_start,
                                  "to", sc_coordinatesObject$refine_stop))

      }
    })
  })


  observeEvent (input$sc_demo_data,{
    isolate({
			 data("singlecell_subset", package="methylscaper", envir = environment())
      sc_seq_data$gch <- singlecell_subset$gch
      sc_seq_data$hcg <- singlecell_subset$hcg
			
			showNotification("Data successfully loaded! Please select Mouse 
			under Choose Organism and select a gene (e.g., Eef1g, Mta2, or Tut1).", type="default", duration=8)
		})
  })
    

  output$sc_seqPlot <- renderPlot({
    obj <- sc_orderObject
    if (sum(obj$toClust) == 0) {}
    else drawPlot(obj,isolate(sc_coordinatesObject))
  }, height = function() {
      session$clientData$output_sc_seqPlot_width
    })


  output$sc_plot_down <- downloadHandler(
    filename = function(){
      if (input$sc_plot_filetype == "PNG") return(paste0("methylscaper_",outname_rds$usename,".png"))
      if (input$sc_plot_filetype == "PDF") return(paste0("methylscaper_",outname_rds$usename,".pdf"))
    },
    content = function(file){
      if (input$sc_plot_filetype == "PNG") png(file)
      if (input$sc_plot_filetype == "PDF") pdf(file)

      drawPlot(sc_orderObject, sc_coordinatesObject, 
                  drawLines = FALSE, plotFast = FALSE)
      dev.off()
    }
  )

  output$sc_log_down <- downloadHandler(
    filename = function(){
      paste0("methylscaper_log_",outname_rds$usename,".txt")
    },
    content = function(file){
      fileConn <- file(file)
      writeLines(actionsLog$log, fileConn)
      close(fileConn)
    }
  )

  output$sc_info <- renderText({
    paste0("Refinement selection: ", sc_coordinatesObject$refine_start, 
                                " ", sc_coordinatesObject$refine_stop, "\n",
           "Weighting selection: ", sc_coordinatesObject$weight_start,
                                " ", sc_coordinatesObject$weight_stop)
  })

	observe({
	    if (sum(sc_orderObject$toClust) == 0) {
	      shinyjs::disable("sc_proportion_hist_download")
			  shinyjs::disable("sc_proportion_data_download")
	      shinyjs::disable("sc_percentC_plot_download")
			  shinyjs::disable("sc_percentC_data_download")
	      shinyjs::disable("sc_avg_c_plot_download")
			  shinyjs::disable("sc_avg_c_data_download")
	      shinyjs::disable("sc_plot_down")
			  shinyjs::disable("sc_log_down")
	    } else {
	      shinyjs::enable("sc_proportion_hist_download")
			  shinyjs::enable("sc_proportion_data_download")
	      shinyjs::enable("sc_percentC_plot_download")
			  shinyjs::enable("sc_percentC_data_download")
	      shinyjs::enable("sc_avg_c_plot_download")
			  shinyjs::enable("sc_avg_c_data_download")
	      shinyjs::enable("sc_plot_down")
			  shinyjs::enable("sc_log_down")}
	  })
		
  output$sc_proportion_color_histogram <- renderPlot({
    obj <- sc_orderObject
    if (sum(obj$toClust) == 0) {}
    else {par(mar=c(5,4,2,2))
			methyl_proportion(obj, makePlot = TRUE,   
            type = input$sc_proportion_choice, main="", xlab="Proportion methylation within cells")
					}
  })

  output$sc_proportion_hist_download <- downloadHandler(
    filename = function(){
        return(paste0("prop_cell_methylated_", tolower(input$sc_proportion_choice), "_", outname_rds$usename, ".pdf"))
    },
    content = function(file){
                  pdf(file)
                  methyl_proportion(sc_orderObject, makePlot = TRUE,
                                   type = input$sc_proportion_choice, main="Methylated Bases Per Cell")
      dev.off()
    }
  )
  output$sc_proportion_data_download <- downloadHandler(
    filename = function(){
        return(paste0("prop_cell_methylated_", tolower(input$sc_proportion_choice), "_", outname_rds$usename, ".csv"))
    },
    content = function(file){
        dat <-  methyl_proportion(sc_orderObject, makePlot = FALSE,
                               type = input$sc_proportion_choice, main="")
      write.csv(dat, file = file)
    }
  )

  output$sc_percent_C <- renderPlot({
    if (sum(sc_orderObject$toClust) == 0) {}
    else {par(mar=c(5,4,2,2))
			methyl_percent_sites(sc_orderObject, makePlot=TRUE)
		}
  })

  output$sc_percentC_plot_download <- downloadHandler(
    filename = function(){
        return(paste0("prcnt_bases_methylated_", outname_rds$usename, ".pdf"))
    },
    content = function(file){
        pdf(file)
        methyl_percent_sites(sc_orderObject, makePlot = TRUE)
      dev.off()
    }
  )

  output$sc_percentC_data_download <- downloadHandler(
    filename = function(){
      return(paste0("prcnt_bases_methylated_", outname_rds$usename, ".txt"))
    },
    content = function(file){
      dat <-  methyl_percent_sites(sc_orderObject, makePlot = FALSE)
      capture.output(dat, file = file)
    }
  )

  output$sc_avg_c <- renderPlot({
    obj <- sc_orderObject
    if (sum(obj$toClust) == 0){}
    else {
			par(mar=c(5,4,2,2))
			methyl_average_status(obj, makePlot=TRUE, window_length=input$sc_window_choice)
		}
  })
	
  output$sc_avg_c_data_download <- downloadHandler(
    filename = function(){
        return(paste0("avg_prcnt_bases_methylated_", outname_rds$usename,".txt"))
    },
    content = function(file){
      dat <-  methyl_average_status(sc_orderObject, makePlot = FALSE, window_length=input$sc_window_choice)
      capture.output(dat, file = file)
    }
  )
	
  output$sc_avg_c_plot_download <- downloadHandler(
    filename = function(){
            return(paste0("avg_prcnt_bases_methylated_", outname_rds$usename,".pdf"))
    },
    content = function(file){
        pdf(file)
        methyl_average_status(sc_orderObject, makePlot = TRUE, window_length=input$sc_window_choice)
      dev.off()
    }
  )
	


  ########################
  # Single-molecule data #
  ########################
  sm_input_data <- reactiveValues(gch = NULL, hcg = NULL)
  sm_raw_data <- reactiveValues(gch = NULL, hcg = NULL)


  observe({
  		if (input$seriate_sm == "Preprocessing" & input$big_tab == "Single-molecule") {
				showNotification("Please provide reference and FASTA files to begin",
                                    type="message", duration=4)
      }
		 if (input$seriate_sm == "Visualization" & input$big_tab == "Single-molecule") {
			showNotification("Provide select the RDS file to begin",
                                  type="message", duration=4)
			}
	})
	
  # alignment handling
  observeEvent(input$run_align, {
    validate(
			need(input$ref_file$datapath, "Please provide the reference .fasta file."),
			need(input$fasta_file$datapath, "Please provide the reads .fasta file.")
			)
			
		 ref <- tryCatch(read.fasta(input$ref_file$datapath), 
		 error=function(cond) {message(paste("Please check the format of your .fasta file"))
		 											# Choose a return value in case of error
                            return(NA)
                        		})
      if (!is.list(ref)) {
              showNotification("Please check the format of your reference .fasta file",
                         type="error", duration=4)
      }      

			
    fasta <- tryCatch(read.fasta(input$fasta_file$datapath),
                       error=function(cond) {
                                   message(paste("Please check the format of your .fasta file"))
                                   # Choose a return value in case of error
                                   return(NA)
                               })
     if (!is.list(fasta)) {
             showNotification("Please check the format of your reads .fasta file",
                        type="error", duration=4)
			}       
     
     if (is.list(ref) & is.list(fasta)) {
             
			if (length(ref)==1){ref <- ref[[1]]}

    progress <- Progress$new()
    progress$set(message = "Beginning alignment", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      progress$set(value = value, message = message, detail = detail)}

    align_out <- tryCatch(runAlign(ref, fasta, updateProgress = updateProgress,
                          log_file = input$processing_log_name),
                               error=function(cond) {
                                           message(paste("No good alignments were found."))
                                           # Choose a return value in case of error
                                           return(NA)
     })

    if (!is.list(align_out)) {
           showNotification("No good alignments were found.",
                      type="error", duration=4)
     }  
    if(is.list(align_out)) {
     sm_raw_data$gch <- align_out$gch
     sm_raw_data$hcg <- align_out$hcg
     sm_raw_data$log_vector  <- align_out$logs
	  }
   read_name <- tools::file_path_sans_ext(input$fasta_file$name)
   ref_name <- tools::file_path_sans_ext(input$ref_file$name)
   outname$usename <- paste0(read_name, "_", ref_name)
  }

 })

	observe({
	    if (is.null(sm_raw_data$hcg)){
				shinyjs::disable("sm_preprocessing_down") 
			  shinyjs::disable("processing_log")
	   } else {
	      shinyjs::enable("sm_preprocessing_down")
			  shinyjs::enable("processing_log")}
	  })
		
		
output$sm_preprocessing_down <- downloadHandler(
  filename = function(){
            paste0(outname$usename,".rds")
    },
    content = function(file){
      saveRDS(list(gch = sm_raw_data$gch, hcg = sm_raw_data$hcg), file = file)
    }
  )

  output$processing_log <- downloadHandler(
    filename = function(){
        paste0("Preprocessing_log_",outname$usename,".txt")
      },
      content = function(file){
        writeLines(sm_raw_data$log_vector, con=file)
      }
    )
    
  observe({if (!is.null(input$sm_rds_file))
  {
    temp <- readRDS(file = input$sm_rds_file$datapath)
    temp_gch <- temp$gch
    temp_hcg <- temp$hcg
		outname_rds$usename <- tools::file_path_sans_ext(input$sm_rds_file$name)
    if (all(rownames(temp_hcg) == temp_hcg[,1])) temp_hcg <- temp_hcg[,-1]
    if (all(rownames(temp_gch) == temp_gch[,1])) temp_gch <- temp_gch[,-1]
		
    if (nrow(temp_gch) == nrow(temp_hcg))
    {
      sm_coordinatesObject$refine_start <- 0
      sm_coordinatesObject$refine_stop <- 0
      sm_coordinatesObject$weight_start <- 0
      sm_coordinatesObject$weight_stop <- 0
      sm_input_data$gch <- temp_gch
      sm_input_data$hcg <- temp_hcg
      sm_input_data$datatype <- "sm"
      isolate({
        actionsLog$log <- c(actionsLog$log, paste("Beginning 
            single-molecule data analysis"))
        actionsLog$log <- c(actionsLog$log, paste("Loading data:",
            input$sm_rds_file$name))
      })
    }

  }})

  # this object keeps track of the coordinates for refinement and weighting
  sm_coordinatesObject <- reactiveValues(refine_start = 0, refine_stop = 0,
                                      weight_start = 0, weight_stop = 0, 
                                      weight_color = "red")
  # now construct the sm_orderObject
  sm_orderObject <- reactiveValues(toClust = 0, order1 = 0)
  observe({ if (!is.null(sm_input_data$gch) & !is.null(sm_input_data$hcg))
  {
   
    progress <- Progress$new()
    progress$set(message = "Beginning seriation", value = 0)
    on.exit(progress$close())

    updateProgress <- function(value = NULL, message = NULL, detail = NULL) {
      progress$set(value = value, message = message, detail = detail)}

    tempObj <- buildOrderObjectShiny(sm_input_data, 
                        input$sm_ser_method, sm_coordinatesObject, updateProgress)
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
    validate(need(sm_input_data$gch, "Please provide input data"))
    validate(need(sm_input_data$hcg, "Please provide input data"))  
    n <- nrow(sm_input_data$gch)
    m <- ncol(sm_input_data$hcg)
    processed_brush <- handleBrushCoordinates(input$sm_plot_brush, n, m)

    if (isolate(input$sm_brush_choice) == "Weighting")
    {
      sm_coordinatesObject$refine_start <- 0
      sm_coordinatesObject$refine_stop <- 0
      sm_coordinatesObject$weight_start <- processed_brush$first_col
      sm_coordinatesObject$weight_stop <- processed_brush$last_col
      sm_coordinatesObject$weight_color <- processed_brush$weight_color
      isolate({
        actionsLog$log <- c(actionsLog$log,
                            paste("Weighting", 
                                    processed_brush$weight_color, "columns",
                                  processed_brush$first_col, "to",
                                  processed_brush$last_col))
      })
    }
    if (isolate(input$sm_brush_choice) == "Refinement")
    {
      s <- processed_brush$first_row
      f <- processed_brush$last_row
      if (s != f)
      {
        sm_coordinatesObject$refine_start <- s
        sm_coordinatesObject$refine_stop <- f
        sm_orderObject$order1 <- refineOrderShiny(isolate(sm_orderObject),
                                    refine_method = isolate(input$sm_refine_method),
                                    sm_coordinatesObject)
        isolate({
          actionsLog$log <- c(actionsLog$log,
                              paste("Refining rows",
                                    processed_brush$first_row, "to",
                                    processed_brush$last_row))
          actionsLog$log <- c(actionsLog$log,
                              paste("Applying refinement with", 
                                      input$sm_refine_method))
        })
      }
    }


  })

  observeEvent( input$sm_force_reverse, {
    isolate({
      if (sm_coordinatesObject$refine_start == sm_coordinatesObject$refine_stop)
      {
        sm_orderObject$order1 <- rev(sm_orderObject$order1)
        actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to", 
                                nrow(sm_input_data$gch)))
      }
      else
      {
        sm_orderObject$order1[sm_coordinatesObject$refine_start : sm_coordinatesObject$refine_stop] <-
          sm_orderObject$order1[sm_coordinatesObject$refine_stop : sm_coordinatesObject$refine_start]
        actionsLog$log <- c(actionsLog$log,
                            paste("Reversing rows", sm_coordinatesObject$refine_start,
                                  "to", sm_coordinatesObject$refine_stop))

      }
    })
  })

  observeEvent(input$sm_demo_data, {
    isolate({
			data("singlemolecule_example", package="methylscaper", envir = environment())
      sm_input_data$gch <- singlemolecule_example$gch
      sm_input_data$hcg <- singlemolecule_example$hcg
    })
  })


  output$sm_seqPlot <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0) {}
    else drawPlot(obj,isolate(sm_coordinatesObject))
  }, height = function() {
      session$clientData$output_sm_seqPlot_width
    })

		

  output$sm_plot_down <- downloadHandler(
    filename = function(){
      if (input$sm_filetype == "PNG") return(paste0("methylscaper_",outname_rds$usename,".png"))
      if (input$sm_filetype == "PDF") return(paste0("methylscaper_",outname_rds$usename,".pdf"))
    },
    content = function(file){
      if (input$sm_filetype == "PNG") png(file)
      if (input$sm_filetype == "PDF") pdf(file)

      drawPlot(sm_orderObject, sm_coordinatesObject, 
                  drawLines = FALSE, plotFast = FALSE)
      dev.off()
    }
  )

  output$sm_log_down <- downloadHandler(
    filename = function(){
      paste0("methylscaper_log_",outname_rds$usename,".txt")
    },
    content = function(file){
      fileConn <- file(file)
      writeLines(actionsLog$log, fileConn)
      close(fileConn)
    }
  )

  output$sm_info <- renderText({
    paste0("Refinement selection: ", sm_coordinatesObject$refine_start, " ",
                 sm_coordinatesObject$refine_stop, "\n",
           "Weighting selection: ", sm_coordinatesObject$weight_start, " ",
                    sm_coordinatesObject$weight_stop)
  })


	observe({
	    if (sum(sm_orderObject$toClust) == 0) {
	      shinyjs::disable("sm_proportion_hist_download")
			  shinyjs::disable("sm_proportion_data_download")
	      shinyjs::disable("sm_percentC_plot_download")
			  shinyjs::disable("sm_percentC_data_download")
	      shinyjs::disable("sm_avg_c_plot_download")
			  shinyjs::disable("sm_avg_c_data_download")
				shinyjs::disable("sm_plot_down") 
			  shinyjs::disable("sm_log_down")
	    } else {
	      shinyjs::enable("sm_proportion_hist_download")
			  shinyjs::enable("sm_proportion_data_download")
	      shinyjs::enable("sm_percentC_plot_download")
			  shinyjs::enable("sm_percentC_data_download")
	      shinyjs::enable("sm_avg_c_plot_download")
			  shinyjs::enable("sm_avg_c_data_download")
	      shinyjs::enable("sm_plot_down")
			  shinyjs::enable("sm_log_down")}
	  })

		
  output$sm_proportion_color_histogram <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0) {}
    else {par(mar=c(5,4,2,2))
			methyl_proportion(obj, makePlot = TRUE, 
                type = input$sm_proportion_choice, main="",xlab="Proportion methylation within molecules")
							}
  })

  output$sm_proportion_hist_download <- downloadHandler(
    filename = function(){
       if (input$sm_proportion_choice == "Accessibility Methylation") {
               whichMeth <- "acc"
       } else {whichMeth <- "met"}
       return(paste0("prop_molecule_methylated_", whichMeth, "_", outname_rds$usename, ".pdf"))
    },
    content = function(file){
       pdf(file)
       methyl_proportion(sm_orderObject, makePlot = TRUE,
                       type = input$sm_proportion_choice, main="Methylated Bases Per Molecule")
      dev.off()
    }
  )
  output$sm_proportion_data_download <- downloadHandler(
    filename = function(){
      if (input$sm_proportion_choice == "Accessibility Methylation") {
              whichMeth <- "acc"
      } else {whichMeth <- "met"}
      return(paste0("prop_molecule_methylated_", whichMeth, "_", outname_rds$usename, ".csv"))
    },
    content = function(file){
      dat <-  methyl_proportion(sm_orderObject, makePlot = FALSE,
                               type = input$sm_proportion_choice, main="Methylated Basepairs Per Molecule")
      write.csv(dat, file = file)
    }
  )

  output$sm_percent_C <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0){}
    else {par(mar=c(5,4,2,2))
			methyl_percent_sites(obj, makePlot=TRUE)
		}
  })

  output$sm_percentC_plot_download <- downloadHandler(
    filename = function(){
            return(paste0("prcnt_bases_methylated_", outname_rds$usename, ".pdf"))
    },
    content = function(file){
        pdf(file)
        methyl_percent_sites(sm_orderObject, makePlot = TRUE)
      dev.off()
    }
  )

  output$sm_percentC_data_download <- downloadHandler(
    filename = function(){
        return(paste0("prcnt_bases_methylated_", outname_rds$usename, ".txt"))
    },
    content = function(file){
      dat <-  methyl_percent_sites(sm_orderObject, makePlot = FALSE)
      capture.output(dat, file = file)
    }
  )

  output$sm_avg_c <- renderPlot({
    obj <- sm_orderObject
    if (sum(obj$toClust) == 0){}
    else {par(mar=c(5,4,2,2))
			methyl_average_status(obj, makePlot=TRUE, window_length=input$sm_window_choice)
		}
  })
	
  output$sm_avg_c_data_download <- downloadHandler(
    filename = function(){
        return(paste0("avg_prcnt_bases_methylated_", outname_rds$usename,".txt"))
    },
    content = function(file){
      dat <-  methyl_average_status(sm_orderObject, makePlot = FALSE, window_length=input$sm_window_choice)
      capture.output(dat, file = file)
    }
  )
	
  output$sm_avg_c_plot_download <- downloadHandler(
    filename = function(){
            return(paste0("avg_prcnt_bases_methylated_", outname_rds$usename,".pdf"))
    },
    content = function(file){
        pdf(file)
        methyl_average_status(sm_orderObject, makePlot = TRUE, window_length=input$sm_window_choice)
      dev.off()
    }
  )
	

}
