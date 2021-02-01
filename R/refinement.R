#' Refinement step
#' 
#' Reorders a subset of the methylation data. 
#' 
#' @param orderObject An object of class \code{orderObject}, generated 
#'  with the \code{initialOrder} function.
#' @param refineStart The index of the first sample (row) used 
#'  in the refinement.
#' @param refineEnd The index of the last sample used in the refinement.
#' @param Method The seriation method used to perform the refinement.
#' 
#' @return The new complete ordering with the refinement applied.
#' @export
#'
#' @examples 
#'  
#' data(day7)
#' 
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")
#' # reorder first 50 cells/molecules (rows)
#' orderObj$order1 <- refineFunction(orderObj, 1, 50) 

refineFunction <- function(orderObject, refineStart, refineEnd, 
                                    Method="HC_average") {
  
  toClust <- orderObject$toClust
  order1 <- orderObject$order1
  toRefine.order <- order1[refineStart:refineEnd]
  
  
  toRefine.clust <- toClust[toRefine.order,]
  
  if (Method=="PCA") {
    col.centered <- apply(toRefine.clust, 2, function(x) x - mean(x))
    try1 <- svd(col.centered, nu = 1, nv = 0)
    order.new <- order(try1$u[,1])
  } else { # Methods available for refining are: ARSA, HC_complete, 
      ## HC_average, HC_ward.
    if (is.null(orderObject$distMat)) distMat <- dist(toRefine.clust,
                                                    method = "euclidean")
    else {
        distMat <- as.dist(as.matrix(orderObject$distMat)[toRefine.order,toRefine.order])
    }
    order.new <- seriation::seriate(distMat, method=Method, verbose=FALSE)
    order.new <- seriation::get_order(order.new)
  }
  
  
  # New order:
  order.new <- order1[seq(refineStart,refineEnd)][order.new]
  order.final <- order1
  order.final[refineStart:refineEnd] <- order.new
  
  return(order.final)
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
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")
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
