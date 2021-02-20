buildOrderObjectShiny <- function(gch, hcg, method,
                        coordinatesObject, updateProgress)
{
    if (coordinatesObject$weight.start == 0 | coordinatesObject$weight.stop == 0) {
        orderObject <- initialOrder(gch, hcg, Method = method)
    } else orderObject <- initialOrder(gch, hcg, Method = method,
                    weightStart = coordinatesObject$weight.start,
                    weightEnd = coordinatesObject$weight.stop,
                    weightFeature = coordinatesObject$weight.color,
                    updateProgress = updateProgress)
return(orderObject)
}

refineOrderShiny <- function(orderObject, refine.method, coordinatesObject)
{
    refineFunction(orderObject, coordinatesObject$refine.start, 
                    coordinatesObject$refine.stop, Method=refine.method)
}

drawPlot <- function(orderObject, coordinatesObject, drawLines = TRUE, ...)
{
    plotSequence(orderObject, ...)
    # draw the horizontal lines
    if (coordinatesObject$refine.start != 0 & coordinatesObject$refine.stop != 0) {
        n <- nrow(orderObject$toClust)# convert back to raw coordinates
        ymin <- (((n:1)[coordinatesObject$refine.start] / n * (n - 10)) + 10) / n 
        ymax <- (((n:1)[coordinatesObject$refine.stop] / n * (n - 10)) + 10) / n
        if (drawLines) {
            abline(b = 0, a = ymax, col = "blue", lwd = 2.5)
            abline(b = 0, a = ymin, col = "blue", lwd = 2.5)
        }
    }
    # draw the vertical lines
    if (coordinatesObject$weight.start != 0 & coordinatesObject$weight.stop != 0) {
        m <- ncol(orderObject$toClust) / 2 # convert back to raw coordinates
        xmin <- (coordinatesObject$weight.start / m) * 0.45
        xmax <- (coordinatesObject$weight.stop / m) * 0.45
        if (coordinatesObject$weight.color == "yellow") {
            xmin <- xmin + 0.55
            xmax <- xmax + 0.55
        }
        if (drawLines) {
            abline(v = xmin, col = "green", lwd = 2.5)
            abline(v = xmax, col = "green", lwd = 2.5)
        }
    }
}

handleBrushCoordinates <- function(plot_brush, n, m){
    weight.color <- "red"

    first.row.raw <- round(plot_brush$ymin * n) - 10
    last.row.raw <- round(plot_brush$ymax * n) - 10
    first.row <- round((first.row.raw / (n - 10)) * n)
    last.row <- round((last.row.raw / (n - 10)) * n)

    if (first.row <= 2) first.row <- 1
    if (last.row >= n - 1) last.row <- n

    if (first.row >= n - 1 | last.row <= 2) {
        first.row <- 0
        last.row <- 0
    }

    first.col <- round(plot_brush$xmin, 2)
    last.col <- round(plot_brush$xmax, 2)

    if (first.col <= 0.45) { # red weighting 
        if (last.col >= 0.45) last.col <- 0.45 # force the last column to be in red
        first.col <- first.col / 0.45
        first.col <- round(first.col * m)
        last.col <- last.col / 0.45
        last.col <- round(last.col * m)
    } else if (first.col >= 0.55) { # yellow weighting
        weight.color <- "yellow"
        first.col <- first.col - 0.55
        last.col <- last.col - 0.55

        first.col <- first.col / 0.45
        first.col <- round(first.col * m)
        last.col <- last.col / 0.45
        last.col <- round(last.col * m)
    } else { # in the middle, just set them to 0
        first.col <- 0
        last.col <- 0
    }

    if (first.col <= 2) first.col <- 1
    if (last.col >= (m - 2)) last.col <- m

    return(list(first.row = ifelse(first.row == 0, 0, (n:1)[first.row]),
            last.row = ifelse(last.row == 0, 0, (n:1)[last.row]),
            first.col = first.col, last.col = last.col,
            weight.color = weight.color))
}

