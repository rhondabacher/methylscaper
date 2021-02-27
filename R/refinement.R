#' Refinement
#' 
#' Reorders a subset of the methylation data. 
#' 
#' @param orderObject An object of class \code{orderObject}, generated 
#'  with the \code{initialOrder} function.
#' @param refineStart The index of the first sample (row) used 
#'  in the refinement.
#' @param refineEnd The index of the last sample (row) used in the refinement.
#' @param Method The seriation method used to perform the refinement.
#' 
#' @return The refinement reorders the cells/molecules (rows) between the indicated
#' start and end positions. The function returns the new complete ordering 
#' with the refinement applied.
#' @export
#'
#' @examples 
#'  
#' data(day7)
#' 
#' orderObj <- initialOrder(day7, Method = "PCA")
#' # reordering the first 50 cells/molecules (rows)
#' orderObj$order1 <- refineFunction(orderObj, 1, 50) 

refineFunction <- function(orderObject, refineStart, refineEnd, 
                                    Method="PCA") {
  
  toClust <- orderObject$toClust
  order1 <- orderObject$order1
  
  if (refineEnd > max(order1) | refineStart > max(order1)) {
     message(paste0("Refine parameters are out of bounds. The maximum has been set to ", refineEnd))
     if (refineEnd > max(order1)) refineEnd <- max(order1)
     if (refineStart > max(order1)) refineStart <- max(order1)  
  }
  if (refineEnd < min(order1) | refineStart < min(order1)) {
     message(paste0("Refine parameters are out of bounds. The maximum has been set to ", refineEnd))
     if (refineEnd < min(order1)) refineEnd <- min(order1)
     if (refineStart < min(order1)) refineStart <- min(order1)  
  }
  
 
  
  toRefine_order <- order1[refineStart:refineEnd]
  
  
  toRefine_clust <- toClust[toRefine_order,]
  
  if (Method=="PCA") {
    col_centered <- apply(toRefine_clust, 2, function(x) x - mean(x))
    try1 <- svd(col_centered, nu = 1, nv = 0)
    order_new <- order(try1$u[,1])
  } else { # Methods available for refining are: ARSA, HC_complete, 
      ## HC_average, HC_ward.
    if (is.null(orderObject$distMat)) distMat <- dist(toRefine_clust,
                                                    method = "euclidean")
    else {
        distMat <- as.dist(as.matrix(orderObject$distMat)[toRefine_order,toRefine_order])
    }
    order_new <- seriation::seriate(distMat, method=Method, verbose=FALSE)
    order_new <- seriation::get_order(order_new)
  }
  
  
  # New order:
  order_new <- order1[seq(refineStart,refineEnd)][order_new]
  order_final <- order1
  order_final[refineStart:refineEnd] <- order_new
  
  return(order_final)
}

#' Force reversal of a subset of the ordering
#' 
#' This reverses a subset of the ordering, as determined by the user. 
#' By default, the entire ordering is reversed.
#' 
#' @param orderObject An object of class \code{orderObject}, 
#'  generated with the \code{initialOrder} function.
#' @param reverseStart The first index to be included in the reversal.
#' @param reverseEnd The last index to be included in the reversal.
#' 
#' @return The new complete ordering, with the reversal applied.
#' @export
#' @examples 
#'  
#' data(day7)
#' 
#' orderObj <- initialOrder(day7, Method = "PCA")
#' # reorder first 50 cells/molecules (rows)
#' orderObj$order1 <- refineFunction(orderObj, 1, 50) 
#' orderObj$order1 <- forceReverse(orderObj, 1, 50) 

forceReverse <- function(orderObject, reverseStart = 1, 
                        reverseEnd = length(orderObject$order1))
{
  order1 <- orderObject$order1
  order1[reverseStart:reverseEnd] <- rev(order1[reverseStart:reverseEnd])
  return(order1)
}
