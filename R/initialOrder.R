#' Ordering of methylation data
#'
#' Orders methylation data using a given seriation method. This function can
#' also perform a weighted seriation if the method is set to "PCA".
#'
#' @param input.GCH The GCH input data.
#' @param input.HCG The HCG input data.
#' @param Method Indicates the seriation method to use. The default option
#'  is "PCA", which orders the data using the first principal component. Any 
#' seriation method provided in the \code{seriation} package is valid.
#' @param weightStart Index of the first column used in the weighted seriation.
#' @param weightEnd Index of the last column used in the weighted seriation.
#' @param weightFeature Indicates whether to weight the GCH or HCG data.
#' @param updateProgress A function to handle the progress bar for the 
#' Shiny app. Should not be used when using the function independently.
#'
#' @return An object of class \code{orderObject}, which contains the generated
#'  ordering and the cleaned data matrix.
#' @importFrom seriation seriate get_order
#' @importFrom stats as.dist dist
#' @importFrom Rfast Dist
#' @export
#' @examples 
#'  
#' data(day7)
#' 
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")


initialOrder <- function(input.GCH, input.HCG, Method="PCA", weightStart=NULL,
                        weightEnd=NULL, weightFeature="red", 
                        updateProgress = NULL){

    ## File checks:
    if (nrow(input.HCG) != nrow(input.GCH)) {
        stop("Input files have different numbers of rows.")}

    if (all(rownames(input.GCH) == input.GCH[,1])) input.GCH <- input.GCH[,-1]
    if (all(rownames(input.HCG) == input.HCG[,1])) input.HCG <- input.HCG[,-1]

    if (is.function(updateProgress)) {
        updateProgress(message = "Recoding input data", value = 0.1)}
    recoded <- recode(input.GCH, input.HCG)
    input.GCH <- recoded$input.GCH
    input.HCG <- recoded$input.HCG

    ## Clustering:
    toClust <- cbind(input.GCH, input.HCG)
    weighted = FALSE

    # (Optional) Weighting: Adds a variable indicating the number 
    # of red or yellow patches at specific DNA location
    if (!is.null(weightStart) & !is.null(weightEnd)) {
        if (is.function(updateProgress)) {
            updateProgress(message = "Weighting selected columns", value = 0.2)
        }
        weighted = TRUE
        if (weightFeature == "red") {
            FEATURE = 3
            weightVector <- apply(input.HCG[,weightStart:weightEnd], 
                                    1, function(x) sum(x==FEATURE))
        }
        if (weightFeature == "yellow") {
            FEATURE = -3
            weightVector <- apply(input.GCH[,weightStart:weightEnd], 
                                    1, function(x) sum(x==FEATURE))
        }
        weightVector[weightVector == 0] <- 1 # we dont want to have 0 weights
    }

    if (is.function(updateProgress)) {
        updateProgress(message = paste("Ordering with", Method), value = 0.35)}
    ## PCA should be the default method:
    if (Method=="PCA") {
        if (weighted) {
            w <- weightVector / sum(weightVector)
            w.sqrt <- sqrt(w)
            toClust.weighted <- diag(w.sqrt) %*% toClust

            col.centered <- apply(toClust.weighted, 2, function(x) x - mean(x))
            try1 <- svd(col.centered, nu = 1, nv = 0)
            order1 <- order(try1$u[,1])
        } else {
            col.centered <- apply(toClust, 2, function(x) x - mean(x))
            try1 <- svd(col.centered, nu = 1, nv = 0)
            order1 <- order(try1$u[,1])
        }
    } else {
        if (weighted) {
            w <- weightVector / sum(weightVector)
            w.sqrt <- sqrt(w)
            toClust.weighted <- diag(w.sqrt) %*% toClust

            distMat <- as.dist(Rfast::Dist(toClust.weighted,
                                    method = "euclidean"))
            order1 <- seriation::seriate(distMat, method=Method)
            order1 <- seriation::get_order(order1)
        } else {
            distMat <- as.dist(Rfast::Dist(toClust,method = "euclidean"))
            order1 <- seriation::seriate(distMat, method=Method)
            order1 <- seriation::get_order(order1)

        }
    }
    orderObject <- list(toClust = toClust, order1 = order1)
    if (Method != "PCA") orderObject$distMat <- distMat
    if (weighted) orderObject$weights <- weightVector
    if (is.function(updateProgress)) {
        updateProgress(message = "Done", value = 1)}
    return(orderObject)
}


recode <- function(input.GCH, input.HCG)
{
    input.GCH[input.GCH=="."] <- 99
    input.GCH <- apply(input.GCH, 2, as.numeric)
    input.GCH[input.GCH==2] <- -4
    input.GCH[input.GCH==1] <- -3
    input.GCH[input.GCH==0] <- -2.5
    input.GCH[input.GCH==-1] <- -99
    input.GCH[input.GCH==-2] <- -1
    input.GCH[input.GCH==-99] <- -2
    input.GCH[input.GCH==99] <- 0

    input.HCG[input.HCG=="."] <- 99
    input.HCG <- apply(input.HCG, 2, as.numeric)
    input.HCG[input.HCG==2] <- 4
    input.HCG[input.HCG==1] <- 3
    input.HCG[input.HCG==0] <- 2.5
    input.HCG[input.HCG==-1] <- 2
    input.HCG[input.HCG==-2] <- 1
    input.HCG[input.HCG==99] <- 0

    return(list(input.GCH = input.GCH, input.HCG = input.HCG))

}
