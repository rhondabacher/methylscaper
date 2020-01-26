#' Refinement step
#' 
#' Reorders a subset of the methylation data. 
#' 
#' @param orderObject An object of class \code{orderObject}, generated with the \code{initialOrder} function.
#' @param refineStart The index of the first sample (row) used in the refinement.
#' @param refineEnd The index of the last sample used in the refinement.
#' @param Method The seriation method used to perform the refinement.
#' 
#' @return The new complete ordering with the refinement applied.
#' @export
refineFunction <- function(orderObject, refineStart, refineEnd,  Method="HC_average") {
  
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
    else distMat <- orderObject$distMat[toRefine.order,toRefine.order]
    order.new <- seriation::seriate(distMat, method=Method, verbose=FALSE)
    order.new <- seriation::get_order(order.new)
  }
  
  
  # New order:
  order.new <- order1[refineStart:refineEnd][order.new]
  order.final <- order1
  order.final[refineStart:refineEnd] <- order.new
  
  
  # if (isTRUE(reverse)) {order.new <- rev(order.new)}
  
  # This is because we only want the updated order of the intended rows:
  # if (refineStart > 1 & refineEnd < nrow(toClust)) {
  #   order.final <- c(order1[seq(1,(refineStart-1))], 
  #                    order.new, 
  #                    order1[seq(from=(refineEnd+1), to=nrow(toClust))])
  # } else if (refineStart == 1) {
  #   order.final <- c(order.new, 
  #                    order1[seq((refineEnd+1),nrow(toClust))])
  # }  else if (refineEnd == nrow(toClust)) {
  #   order.final <- c(order1[seq(1,(refineStart-1))],
  #                    order.new)
  # }
  return(order.final)
}

