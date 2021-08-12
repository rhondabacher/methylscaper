#' Ordering the molecules/reads
#'
#' This function performs the weighted seriation procedure
#' described in the methylscaper manuscript if the method is set to "PCA".
#' The data may also be ordered using a given seriation method from the
#' seriation R package. The weighting is done between a designated start and end
#' base pair chosen by the user, and the weight can be done on the 
#' endogenous methylation or the accessibility.
#'
#' @param dataIn A list object containing two elements labelled gch and hcg (already pre-processed.)
#' @param Method Indicates the seriation method to use. The default option
#'  is "PCA", which orders the data using a weighted first principal component approach. Any 
#'  seriation method provided in the \code{seriation} package is also valid input. 
#' @param weightStart Index of the first column used in the weighted seriation.
#' @param weightEnd Index of the last column used in the weighted seriation.
#' @param weightFeature Indicates whether to weight the GCH or HCG data.
#' Valid input to weight the GCH is 'gch', 'acc', or 'yellow'. To weight the HCG, 
#' valid input for this option is 'hcg', 'met', or 'red'.
#' @param updateProgress A function to handle the progress bar for the 
#' Shiny app. Should not be used when using the function independently.
#'
#' @return An object of class \code{orderObject}, which contains the generated
#'  ordering ($order1) and a clean data matrix ($toClust) 
#'  to be passed into the plotting function plotSequence().
#' @importFrom seriation seriate get_order
#' @importFrom stats as.dist dist
#' @importFrom Rfast Dist
#' @export
#' @examples 
#'  
#' data(singlemolecule_example)
#' 
#' orderObj <- initialOrder(singlemolecule_example)


initialOrder <- function(dataIn, Method="PCA", weightStart=NULL,
                        weightEnd=NULL, weightFeature="red", 
                        updateProgress = NULL){

    input_GCH <- dataIn$gch 
    input_HCG <- dataIn$hcg
    
    ## File checks:
		
		if(is.null(updateProgress) & !is.list(dataIn)) {
			stop("No valid sites in designated range. Choose another gene or adjust  
			                      start and end positions with a larger range")}
    if (nrow(input_HCG) != nrow(input_GCH)) {
        stop("Input files have different numbers of rows.")}

    if (all(rownames(input_GCH) == input_GCH[,1])) input_GCH <- input_GCH[,-1]
    if (all(rownames(input_HCG) == input_HCG[,1])) input_HCG <- input_HCG[,-1]
				
    if (is.function(updateProgress)) {
        updateProgress(message = "Recoding input data", value = 0.1)}
    
    recoded <- recode(input_GCH, input_HCG)
    input_GCH <- recoded$input_GCH
    input_HCG <- recoded$input_HCG

    weightFeature <- tolower(weightFeature)
    
    
    if (nrow(dataIn$gch) == 1) {
      toClust <- t(as.matrix(c(input_GCH, input_HCG)))
      order1 <- 1
      orderObject <- list(toClust = toClust, order1 = order1)
    } 
    else {
      toClust <- cbind(input_GCH, input_HCG)
      weighted = FALSE

    # (Optional) Weighting: Adds a variable indicating the number 
    # of red or yellow patches at specific DNA location
    if (!is.null(weightStart) & !is.null(weightEnd)) {
        if (is.function(updateProgress)) {
            updateProgress(message = "Weighting selected columns", value = 0.2)
        }
        weighted = TRUE

		if (weightStart < 1) weightStart <- 1
	    if (weightEnd > ncol(input_HCG)) {
			print("Setting weightEnd to maximum columns")
			weightEnd <- ncol(input_HCG)
		}
        if (weightFeature == "red" | weightFeature == 'hcg' | weightFeature == 'met') {
            FEATURE = 3
            weightVector <- apply(input_HCG[,seq(weightStart,weightEnd)], 
                                    1, function(x) sum(x==FEATURE))

        } else if (weightFeature == "yellow" | weightFeature == 'gch' | weightFeature == 'acc') {
           FEATURE = -3
            weightVector <- apply(input_GCH[,seq(weightStart,weightEnd)], 
                                    1, function(x) sum(x==FEATURE))
        } else {
        	print("weightFeature value is not valid, see ?weightFeature for valid values")
        }
        weightVector[weightVector == 0] <- 1 # we dont want to have 0 weights
    }

    if (is.function(updateProgress)) {
        updateProgress(message = paste("Ordering with", Method), value = 0.35)}
    ## PCA should be the default method:
    if (Method=="PCA") {
        if (weighted) {
            w <- weightVector / sum(weightVector)
            w_sqrt <- sqrt(w)
            toClust_weighted <- diag(w_sqrt) %*% toClust

            col_centered <- apply(toClust_weighted, 2, function(x) x - mean(x))
            try1 <- svd(col_centered, nu = 1, nv = 0)
            order1 <- order(try1$u[,1])
        } else {
            col_centered <- apply(toClust, 2, function(x) x - mean(x))
            try1 <- svd(col_centered, nu = 1, nv = 0)
            order1 <- order(try1$u[,1])
        }
    } else {
        if (weighted) {
            w <- weightVector / sum(weightVector)
            w_sqrt <- sqrt(w)
            toClust_weighted <- diag(w_sqrt) %*% toClust

            distMat <- as.dist(Dist(toClust_weighted,
                                    method = "euclidean"))
            order1 <- seriate(distMat, method=Method)
            order1 <- get_order(order1)
        } else {
            distMat <- as.dist(Dist(toClust,method = "euclidean"))
            order1 <- seriate(distMat, method=Method)
            order1 <- get_order(order1)

        }
    }
      orderObject <- list(toClust = toClust, order1 = order1)
      if (weighted) orderObject$weights <- weightVector
      if (Method != "PCA") orderObject$distMat <- distMat
    }
    
    if (is.function(updateProgress)) {
        updateProgress(message = "Done", value = 1)}
    return(orderObject)
}


recode <- function(input_GCH, input_HCG)
{
    input_GCH[input_GCH=="."] <- 99
    input_GCH <- apply(input_GCH, 2, as.numeric)
    input_GCH[input_GCH==2] <- -4
    input_GCH[input_GCH==1] <- -3
    input_GCH[input_GCH==0] <- -2.5
    input_GCH[input_GCH==-1] <- -99
    input_GCH[input_GCH==-2] <- -1
    input_GCH[input_GCH==-99] <- -2
    input_GCH[input_GCH==99] <- 0

    input_HCG[input_HCG=="."] <- 99
    input_HCG <- apply(input_HCG, 2, as.numeric)
    input_HCG[input_HCG==2] <- 4
    input_HCG[input_HCG==1] <- 3
    input_HCG[input_HCG==0] <- 2.5
    input_HCG[input_HCG==-1] <- 2
    input_HCG[input_HCG==-2] <- 1
    input_HCG[input_HCG==99] <- 0


    # Recode the plot edges as white always:
     bp <- length(input_HCG[1,])

     sites = which(apply(input_HCG, 2, function(x) any(x %in% c(4, 1))))
     firstHCG <- sites
     lastHCG <- rev(sites)[1]

     input_HCG <- apply(input_HCG, 1, function(x) {
       x[seq(1,firstHCG)] <- 0
       x[seq(lastHCG, bp)] <- 0
     return(x)
     })
     input_HCG <- t(input_HCG)

     sites = which(apply(input_GCH, 2, function(x) any(x %in% c(-4, -1))))
     firstGCH <- sites
     lastGCH <- rev(sites)[1]

     input_GCH <- apply(input_GCH, 1, function(x) {
       x[seq(1,firstGCH)] <- 0
       x[seq(lastGCH, bp)] <- 0
       return(x)
     })
     input_GCH <- t(input_GCH)
       
    return(list(input_GCH = data.matrix(input_GCH), input_HCG = data.matrix(input_HCG)))

}
