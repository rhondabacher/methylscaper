current.order <- 1:10
unordered <- TRUE
refine.start.g <- 0
refine.stop.g <- 0
weight.start.g <- 0
weight.stop.g <- 0
weight.color.g <- "red"

buildOrderObjectShiny <- function(gch, hcg, method,
                        weight.start, weight.stop, weight.color)
{
  if (weight.start == 0 | weight.stop == 0)
  {
    orderObject <- initialOrder(gch, hcg, Method = method)
  }
  else orderObject <- initialOrder(gch, hcg, Method = method, 
                                   weightStart = weight.start, weightEnd = weight.stop,
                                   weightFeature = weight.color)
  if (unordered) # only save this initial ordering as global if this is the first time the data is ordered
  { # note, this should also happen when a new weighting is assigned. how to check for this?
    current.order <<- orderObject$order1
    unordered <<- FALSE
  }
  orderObject$order1 <- current.order
  return(orderObject)
  
}

refineOrderShiny <- function(orderObject, refine.method)
{
  if (refine.start.g != 0 & refine.stop.g != 0) # save the global ordering to this new refined ordering
    current.order <<- refineFunction(orderObject, refine.start.g, refine.stop.g, refine.method)
}

makePlot <- function(orderObject, ...)
{
 
  orderObject$order1 <- current.order 
  plotSequence(orderObject, ...)
  if (refine.start.g != 0 & refine.stop.g != 0) # draw the horizontal lines
  {
    n <- nrow(orderObject$toClust)
    ymin <- (((140:1)[refine.start.g] / n * (n - 10)) + 10) / n # convert back to raw coordinates
    ymax <- (((140:1)[refine.stop.g] / n * (n - 10)) + 10) / n
    abline(b = 0, a = ymax, col = "blue", lwd = 2.5)
    abline(b = 0, a = ymin, col = "blue", lwd = 2.5)
    
  }
  if (weight.start.g != 0 & weight.stop.g != 0) # draw the vertical lines
  {
    m <- ncol(orderObject$toClust) / 2 # convert back to raw coordinates
    xmin <- (weight.start.g / m) * 0.45
    xmax <- (weight.stop.g / m) * 0.45
    if (weight.color.g == "yellow")
    {
      xmin <- xmin + 0.55
      xmax <- xmax + 0.55
    }
    abline(v = xmin, col = "green", lwd = 2.5)
    abline(v = xmax, col = "green", lwd = 2.5)
  }
}

handleBrushCoordinates <- function(plot_brush, brush.choice, n, m){
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
  
  if (brush.choice == "Refinement")
  {
    refine.start.g <<- (140:1)[first.row]
    refine.stop.g <<- (140:1)[last.row]
  }
  else if (brush.choice == "Weighting")
  {
    unordered <<- TRUE # set the global data to unordered if the weighting changes
    refine.start.g <<- 0
    refine.stop.g <<- 0
    weight.start.g <<- first.col
    weight.stop.g <<- last.col
    weight.color.g <<- weight.color
  }
  
  
  return(list(first.row = (140:1)[first.row],
              last.row = (140:1)[last.row],
              first.col = first.col,
              last.col = last.col,
              weight.color = weight.color))
  
}