## Upload fields:
## input.GCH: The GCH csv from Alberto.
## input.HCG: The HCG csv from Alberto.

# User optional fields:
## weightStart: Where along the DNA to begin the weight the read ordering.
## weightEnd: Where along the DNA to end weight the read ordering.
## weightFeature: Whether to weight the methylation or accessibility ordering.
## reverse: Force the ordering to be exactly reversed.

# returns an initial (unrefined) ordering and to 'toClust' object
#' @import seriation
#' @importFrom aroma.light wpca
initialOrder <- function(input.GCH, input.HCG, Method="PCA", weightStart=NULL, weightEnd=NULL, 
                         weightFeature="red", reverse=F){
  
  ## File checks:
  if (nrow(input.HCG) != nrow(input.GCH)) {stop("Input files have different numbers of rows.")}
  
  # Currently Alberto has a software that is not published that gets us the CSV files.
  # I think we could move most of that to R, but will still need to have a python script 
  # that takes in REF and SEQ's and does the overlap alignment.
  # After that, the rest could ALL be done in R and the recoding would not be necesary.
  ## To do: Talk to Alberto.
  
  # Recoding for heatmaps plots; might be able to move this:
  input.GCH[input.GCH=="."] <- 99
  input.GCH <- data.matrix(input.GCH)
  input.GCH[input.GCH==2] <- -4
  input.GCH[input.GCH==1] <- -3
  input.GCH[input.GCH==0] <- -2.5
  input.GCH[input.GCH==-1] <- -99
  input.GCH[input.GCH==-2] <- -1
  input.GCH[input.GCH==-99] <- -2
  input.GCH[input.GCH==99] <- 0
  
  input.HCG[input.HCG=="."] <- 99
  input.HCG <- data.matrix(input.HCG)
  input.HCG[input.HCG==2] <- 4
  input.HCG[input.HCG==1] <- 3
  input.HCG[input.HCG==0] <- 2.5
  input.HCG[input.HCG==-1] <- 2
  input.HCG[input.HCG==-2] <- 1
  input.HCG[input.HCG==99] <- 0
  
  ## Clustering:  
  toClust <- cbind(input.GCH, input.HCG)
  weighted = FALSE
  
  # (Optional) Weighting: Adds a variable indicating the number of red or yellow patches at specific DNA location
  if (!is.null(weightStart) & !is.null(weightEnd)) {
    weighted = TRUE
    if (weightFeature == "red") {
      FEATURE = 3
      weightVector <- apply(input.HCG[,weightStart:weightEnd], 1, function(x) sum(x[x==FEATURE]))
    }
    if (weightFeature == "yellow") {
      FEATURE = -3
      weightVector <- apply(input.GCH[,weightStart:weightEnd], 1, function(x) sum(x[x==FEATURE]))
    }
    weightVector[weightVector == 0] <- 1 # we dont want to have 0 weights
    weightVector <- abs(weightVector) # weights should only be positive
    # weightMat <- diag(weightVector) # creates the diagonal matrix of weights
    # toClust <- cbind(weightVector, input.GCH, input.HCG) # ask about this, are we accounting for weights in the seriation correctly?
    # toClust <- weightMat %*% cbind(input.GCH, input.HCG)
  }
  
  
  ## PCA should be the default method:
  if (Method=="PCA") {
    if (weighted)
    {
      try1 <- aroma.light::wpca(x = toClust, w = weightVector, center = TRUE)
      order1 <- order(try1$pc[,1])
    }
    else
    {
      col.centered <- apply(toClust, 2, function(x) x - mean(x))
      try1 <- svd(col.centered, nu = 1, nv = 0)
      order1 <- order(try1$u[,1])
    }
    
  } else{ 
    # How do we incorprate the weights here (i.e., when non-PCA methods are used)?
    distMat <- dist(toClust,method = "euclidean") # put in my faster dist code i made before
    # Allow drop down methods to be: ARSA.
    order1 <- seriation::seriate(distMat, method=Method)
    order1 <- seriation::get_order(order1)
  }
  if (isTRUE(reverse)) {order1 <- rev(order1)}
  orderObject <- list(toClust = toClust, order1 = order1)
  if (Method != "PCA") orderObject$distMat <- distMat
  return(orderObject)
}