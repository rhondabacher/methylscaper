#' Calculates the percentage of methylated cells/molecules per site
#'
#' @param orderObject An object of class \code{orderObject}
#' @param makePlot Logical, indicates whether to generate the percentage plot
#' @param ... Additional parameters used by the \code{plot} function.
#' @return The percent of molecules or cells methylated (endogenous (yellow) 
#'      or accessible) at each site. Output is a list with names "red" 
#'      and "yellow". Red represents the endogenous methylation and 
#'      yellow represents the accessibility. Within each list object is a 
#'      vector of the percent of cells/molecules methylated. The location 
#'      of the site is also represent in the form CXX, where XX is the 
#'      position of the site within the defined region.
#' @importFrom graphics hist lines plot points
#' @export
#' @examples 
#'  
#' data(singlemolecule_example)
#' 
#' orderObj <- initialOrder(singlemolecule_example, Method = "PCA")
#' methyl_percent_sites(orderObj, makePlot = TRUE)

methyl_percent_sites <- function(orderObject, makePlot = TRUE, ...){
    dat <- orderObject$toClust
    
    red_sites <- c(which(dat == 4, arr.ind = TRUE)[,2], which(dat == 1, arr.ind = TRUE)[,2])
    red_sites <- sort(unique(red_sites))
    yellow_sites <- c(which(dat == -4, arr.ind = TRUE)[,2], which(dat == -1, arr.ind = TRUE)[,2])
    yellow_sites <- sort(unique(yellow_sites))
    
    c_red <- vapply(red_sites, function(i) {
        sum(dat[, i] == 4) / nrow(dat)}, numeric(1))
    c_yellow <- vapply(yellow_sites, function(i) {
        sum(dat[, i] == -4) / nrow(dat)}, numeric(1))
    if (makePlot) {
        plot(x = red_sites - ncol(dat)/2, y = c_red,
            col="brown1", pch=19, ylim=c(0,1), 
            xlab = "Region (Base Pair)", ylab="%Methylated",
            bty='n', cex.lab=1.3,...)
        lines(x = red_sites - ncol(dat)/2, y = c_red, col="brown1")
        points(x = yellow_sites, y = c_yellow, col="gold2", pch=19)
        lines(x = yellow_sites, y = c_yellow, col="gold2")
        legend("topright", c("Endogenous", "Accessibility"), pch=16, 
                col=c("brown1", "gold2"), bty='n', title="Methylation Type")
        n_sites <- length(union(red_sites, yellow_sites))
        labs <- union(red_sites, yellow_sites)[seq(1, n_sites, by=n_sites/12)]
    }
    final <- list(c_red, c_yellow)
    names(final) = c("meth", "acc")
    return(final)
}

#' Calculate the proportion of methylated bases for each cell/molecule
#'
#' @param orderObject An object of class \code{orderObject}
#' @param type Indicates which data set to compute proportions for.
#'      This should be 'met' or 'hcg' or 'red' for endogenous methylation; 'acc' or
#'      'gch' or 'yellow' for accessibility. 
#' @param makePlot Indicates whether to plot a histogram of the proportions 
#'  across all cells/molcules.
#' @param ... Additional parameters used by the \code{hist} function.
#'
#' @return The proportion of methylated (endogenous (yellow) 
#'      or accessible) bases for each cell/molecule. Output is vector
#'      with length the numbner of cells/molecules and contains a proportion.     
#' @importFrom graphics hist
#' @export
#' @examples 
#'  
#' data(singlemolecule_example)
#' 
#' orderObj <- initialOrder(singlemolecule_example, Method = "PCA")
#' methyl_proportion(orderObj, makePlot = TRUE)

methyl_proportion <- function(orderObject, type = "yellow", 
                                makePlot=TRUE, ...){
  type <- tolower(type)
  if (type=="gch") type <- "yellow"
  if (type=="accessibility methylation") type <- "yellow"
  if (type=="acc") type <- "yellow"
  color_indicator <- ifelse(type=="yellow", -1, 1)
  Proportion <- apply(orderObject$toClust, 1, function(x){
      sum(x == color_indicator * 3 | x == color_indicator * 4) / (length(x) / 2)
  })
  if (makePlot) {
    opar <- par(lwd=4)
    H = hist(Proportion, plot=FALSE, breaks = 15)
    plot(H, xlim=c(0,1), border=ifelse(type == "yellow", "gold2", "brown1"),
        col="gray75", cex.lab=1.3,
        lwd=2,...)
    par(opar)
  }
  Proportion <- data.frame(methylationProportion = Proportion)
  return(Proportion)
}

#' Calculate the average methylation/accessibility status across all cells/molecules.
#'
#' @param orderObject An object of class \code{orderObject}
#' @param window_length Length of the window to be used to compute a 
#'          moving average. Default is 20.
#' @param makePlot Logical, indicates whether to generate a line plot of 
#'          average status.
#' @param ... Addition parameters used by the \code{plot} function.
#'
#' @return The proportion of methylated bases for each cell/molecule 
#'      within a defined moving window. Output is a list with elements 
#'      "meth_avg" and "acc_avg", indicating endogenous
#'      or accessible methylation respectively.
#' @importFrom stats filter
#' @importFrom graphics legend
#' @export
#' @examples 
#'  
#' data(singlemolecule_example)
#' 
#' orderObj <- initialOrder(singlemolecule_example, Method = "PCA")
#' methyl_average_status(orderObj, makePlot = TRUE)
#' 
methyl_average_status <- function(orderObject, window_length = 20, 
                                makePlot = TRUE, ...)
{
    colLength <- ncol(orderObject$toClust)
    gch_num <- orderObject$toClust[,seq(1,(colLength / 2))]
    hcg_num <- orderObject$toClust[,seq((colLength / 2 + 1),colLength)]

    acc_sum <- colSums(gch_num <= -3)
    meth_sum <- colSums(hcg_num >= 3)
    acc_denom <- colSums(gch_num != 0)
    meth_denom <- colSums(hcg_num != 0)
    acc_avg <- acc_sum / acc_denom
    meth_avg <- meth_sum / meth_denom

    width <- window_length
    moving_acc_avg <- filter(x = acc_avg, filter = rep(1, width)) / width
    moving_meth_avg <- filter(x = meth_avg, filter = rep(1, width)) / width

    if (makePlot) {
        plot(x=seq(1,length(moving_acc_avg)), y=moving_acc_avg, 
						type = "l", col="gold2", lwd=2, bty='n',
            xlab="Position along the genomic location", 
            ylab="Averaged methylation status", 
						cex.lab=1.3, ylim = c(0,1),...)
        lines(moving_meth_avg, col="brown1", lwd=2)
        legend("topright", c("Endogenous", "Accessibility"), lwd=2, 
                col=c("brown1", "gold2"), bty='n', title="Methylation Type")
    }
    return(list(meth_avg = moving_meth_avg, acc_avg = moving_acc_avg))
}
