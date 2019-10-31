makePlot <- function(gch, hcg, method, refine.start, refine.stop, refine.method, 
                     weight.start, weight.stop, weight.color, ...)
{
  if (weight.start == 0 | weight.stop == 0)
  {
    orderObject <- initialOrder(gch, hcg, Method = method)
  }
  else orderObject <- initialOrder(gch, hcg, Method = method, 
                              weightStart = weight.start, weightEnd = weight.stop,
                              weightFeature = weight.color)
  
  if (refine.start != 0 & refine.stop != 0) {
    new.order <- refineFunction(orderObject, refineStart = refine.start, refineEnd = refine.stop, 
                                Method = refine.method)
    orderObject$order1 <- new.order
  }
  plotSequence(orderObject, ...)
  if (refine.start != 0 & refine.stop != 0) # draw the horizontal lines
  {
    n <- nrow(gch)
    ymin <- (((140:1)[refine.start] / n * (n - 10)) + 10) / n # convert back to raw coordinates
    ymax <- (((140:1)[refine.stop] / n * (n - 10)) + 10) / n
    abline(b = 0, a = ymax, col = "blue", lwd = 2.5)
    abline(b = 0, a = ymin, col = "blue", lwd = 2.5)
    
  }
  if (weight.start != 0 & weight.stop != 0) # draw the vertical lines
  {
    m <- ncol(gch) # convert back to raw coordinates
    xmin <- (weight.start / m) * 0.45
    xmax <- (weight.stop / m) * 0.45
    if (weight.color == "yellow")
    {
      xmin <- xmin + 0.55
      xmax <- xmax + 0.55
    }
    abline(v = xmin, col = "green", lwd = 2.5)
    abline(v = xmax, col = "green", lwd = 2.5)
  }
}

handleBrushCoordinates <- function(plot_brush, n, m){
  weight.color <- "red"
  
  first.row.raw <- round(plot_brush$ymin * n) - 10
  last.row.raw <- round(plot_brush$ymax * n) - 10
  first.row <- round((first.row.raw / (n - 10)) * n)
  last.row <- round((last.row.raw / (n - 10)) * n)
  
  if(first.row <= 2) first.row <- 1
  if (last.row >= n - 1) last.row <- n
  
  first.col <- round(plot_brush$xmin, 2)
  last.col <- round(plot_brush$xmax, 2)
  
  if (first.col <= 0.45) # red weighting
  {
    if (last.col >= 0.45) last.col <- 0.45 # force the last column to be in the red
    first.col <- first.col / 0.45
    first.col <- round(first.col * m)
    last.col <- last.col / 0.45
    last.col <- round(last.col * m)
    
  }
  else if (first.col >= 0.55) # yellow weighting
  {
    weight.color <- "yellow"
    first.col <- first.col - 0.55
    last.col <- last.col - 0.55
    
    first.col <- first.col / 0.45
    first.col <- round(first.col * m)
    last.col <- last.col / 0.45
    last.col <- round(last.col * m)
  }
  else # in the middle, just set them to 0
  {
    first.col <- 0
    last.col <- 0
  }
  
  if (first.col <= 2) first.col <- 1
  if (last.col >= (m - 2)) last.col <- m
  
  
  return(list(first.row = (140:1)[first.row], # we have to use this indexing, because the coordinates are swapped in the plot
              last.row = (140:1)[last.row],  #
              first.col = first.col,
              last.col = last.col,
              weight.color = weight.color))
  
}