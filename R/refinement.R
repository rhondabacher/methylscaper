#' Refinement step
#' 
#' Reorders a subset of the methylation data. 
#' 
#' @param orderObject An object of class \code{orderObject}, generated with the \code{initialOrder} function.
#' @param refineStart The index of the first sample (row) used in the refinement.
#' @param refineEnd The index of the last sample used in the refinement.
#' @param reverse Logical, indicates whether to reverse the refined ordering.
#' @param Method The seriation method used to perform the refinement.
#' 
#' @return The new complete ordering with the refinement applied.
#' @export
refineFunction <- function(orderObject, refineStart, refineEnd,  reverse = FALSE, Method="HC_average") {
  
  toClust <- orderObject$toClust
  order1 <- orderObject$order1
  toRefine.order <- order1[refineStart:refineEnd]
  
  
  toRefine.clust <- toClust[toRefine.order,]
  
  if (Method=="PCA") {
    col.centered <- apply(toRefine.clust, 2, function(x) x - mean(x))
    try1 <- svd(col.centered, nu = 1, nv = 0)
    order.new <- order(try1$u[,1])
  } else { # Methods available for refining are: ARSA, HC_complete, HC_average, HC_ward.
    ## Shouldn't need to recalculate the distance each time!!
    if (is.null(orderObject$distMat)) distMat <- dist(toRefine.clust,method = "euclidean")
    else distMat <- as.dist(as.matrix(orderObject$distMat)[toRefine.order,toRefine.order])
    order.new <- seriation::seriate(distMat, method=Method, verbose=FALSE)
    order.new <- seriation::get_order(order.new)
  }
  
  
  # New order:
  order.new <- order1[refineStart:refineEnd][order.new]
  order.final <- order1
  order.final[refineStart:refineEnd] <- order.new
  
  if (reverse) order.final <- rev(order.final)
  
  return(order.final)
}

