buildOrderObjectShiny <- function(dataIn, method,
                        coordinatesObject, updateProgress)
{
    if (coordinatesObject$weight_start == 0 | coordinatesObject$weight_stop == 0) {
        orderObject <- initialOrder(dataIn, Method = method)
    } else orderObject <- initialOrder(dataIn, Method = method,
                    weightStart = coordinatesObject$weight_start,
                    weightEnd = coordinatesObject$weight_stop,
                    weightFeature = coordinatesObject$weight_color,
                    updateProgress = updateProgress)
return(orderObject)
}

refineOrderShiny <- function(orderObject, refine_method, coordinatesObject)
{
    refineFunction(orderObject, coordinatesObject$refine_start, 
                    coordinatesObject$refine_stop, Method=refine_method)
}

drawPlot <- function(orderObject, coordinatesObject, blankWidth = 75, drawLines = TRUE, ...)
{
    plotSequence(orderObject, ...)
    # draw the horizontal lines
    if (coordinatesObject$refine_start != 0 & coordinatesObject$refine_stop != 0) {
        n <- nrow(orderObject$toClust)# convert back to raw coordinates
        ymin <- (((n:1)[coordinatesObject$refine_start] / n * (n - 10)) + 10) / n 
        ymax <- (((n:1)[coordinatesObject$refine_stop] / n * (n - 10)) + 10) / n
        if (drawLines) {
            abline(b = 0, a = ymax, col = "blue", lwd = 2.5)
            abline(b = 0, a = ymin, col = "blue", lwd = 2.5)
        }
    }
    # draw the vertical lines
    if (coordinatesObject$weight_start != 0 & coordinatesObject$weight_stop != 0) {
		firstm <- (ncol(orderObject$toClust) + blankWidth)
        xmin <- coordinatesObject$weight_start / firstm
        xmax <- coordinatesObject$weight_stop / firstm
        if (coordinatesObject$weight_color == "yellow") {
			secondm <- (blankWidth+ncol(orderObject$toClust)/2) / firstm
            xmin <- xmin + secondm
            xmax <- xmax + secondm
        }
        if (drawLines) {
            abline(v = xmin, col = "green", lwd = 2.5)
            abline(v = xmax, col = "green", lwd = 2.5)
        }
    }
}

handleBrushCoordinates <- function(plot_brush, n, m, blankWidth=75){
    weight_color <- "red"
    first_row_raw <- round(plot_brush$ymin * n) - 10
    last_row_raw <- round(plot_brush$ymax * n) - 10
    first_row <- round((first_row_raw / (n - 10)) * n)
    last_row <- round((last_row_raw / (n - 10)) * n)

    if (first_row <= 2) first_row <- 1
    if (last_row >= n - 1) last_row <- n

    if (first_row >= n - 1 | last_row <= 2) {
        first_row <- 0
        last_row <- 0
    }

    first_col <- round(plot_brush$xmin, 4)
    last_col <- round(plot_brush$xmax, 4)
	firstm <- (m) / (m*2 + blankWidth)
	secondm <- (blankWidth+m) / (m*2 + blankWidth)
    if (first_col <= firstm) { # red weighting 
        if (last_col >= firstm) last_col <- firstm # force the last column to be in red
        first_col <- first_col / firstm
        first_col <- round(first_col * m)
        last_col <- last_col / firstm
        last_col <- round(last_col * m)
    } else if (first_col >= secondm) { # yellow weighting
        weight_color <- "yellow"
        first_col <- first_col - secondm
        last_col <- last_col - secondm

        first_col <- first_col / firstm
        first_col <- round(first_col * m)
        last_col <- last_col / firstm
        last_col <- round(last_col * m)
    } else { # in the middle, just set them to 0
        first_col <- 0
        last_col <- 0
    }

    if (first_col <= 2) first_col <- 1
    if (last_col >= (m - 2)) last_col <- m

    return(list(first_row = ifelse(first_row == 0, 0, seq(n,1)[first_row]),
            last_row = ifelse(last_row == 0, 0, seq(n,1)[last_row]),
            first_col = first_col, last_col = last_col,
            weight_color = weight_color))
}

