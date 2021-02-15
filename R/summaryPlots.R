#' Calculates the percentage of methylated cells/molecules per site
#'
#' @param orderObject An object of class \code{orderObject}
#' @param makePlot Logical, indicates whether to generate the percentage plot
#' @param ... Additional parameters used by the \code{plot} function.
#' @return The percent of reads or cells methylated (endogenous (yellow) 
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
#' data(day7)
#' 
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")
#' methyl_percent_site(orderObj, plotPercents = TRUE)

methyl_percent_site <- function(orderObject, makePlot = TRUE, ...){
    dat <- orderObject$toClust
    red.sites <- which(dat[1,seq(1,ncol(dat))] == 4 |
                          dat[1,seq(1,ncol(dat))] == 1)

    yellow.sites <- which(dat[1,seq(1,ncol(dat))] == -4 |
                          dat[1,seq(1,ncol(dat))] == -1)
    c.red <- vapply(red.sites, function(i) {
        sum(dat[, i] == 4) / nrow(dat)}, numeric(1))
    c.yellow <- vapply(yellow.sites, function(i) {
        sum(dat[, i] == -4) / nrow(dat)}, numeric(1))
    if (makePlot) {
        plot(x = red.sites - ncol(dat)/2, y = c.red,
            col="brown1", pch=19, ylim=c(0,1), 
            xlab = "Region (Base Pair)", ylab="%Methylated",
            bty='n', cex.lab=1.3, xaxt='n', yaxt='n',...)
        axis(side = 1, lwd = 2, cex.axis=1.2)
        axis(side = 2, lwd = 2, cex.axis=1.2)
        lines(x = red.sites - ncol(dat)/2, y = c.red, col="brown1")
        points(x = yellow.sites, y = c.yellow, col="gold2", pch=19)
        lines(x = yellow.sites, y = c.yellow, col="gold2")
        legend("topright", c("Endogenous", "Accessibility"), pch=16, 
                col=c("brown1", "gold2"), bty='n', title="Methylation Type")
        n.sites <- length(union(red.sites, yellow.sites))
        labs <- union(red.sites, yellow.sites)[seq(1, n.sites, by=n.sites/12)]
    }
    final <- list(c.red, c.yellow)
    names(final) = c("meth", "acc")
    return(final)
}

#' Calculate the proportion of methylated bases for each cell/molecule
#'
#' @param orderObject An object of class \code{orderObject}
#' @param color Indicates which data set to compute proportions for.
#'      This should be 'red' or 'hcg' for endogenous methylation; 'yellow' or
#'      'gch' for accessibility. 
#' @param makePlot Indicates whether to plot a histogram of the proportions 
#'  across all reads.
#' @param ... Additional parameters used by the \code{hist} function.
#'
#' @return The proportion of methylated (endogenous (yellow) 
#'      or accessible) bases for each cell/molecule. Output is vector
#'      with length the numbner of cells/molecules and contains a proportion.     
#' @importFrom graphics hist
#' @export
#' @examples 
#'  
#' data(day7)
#' 
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")
#' methyl_proportion_cell(orderObj, plotHistogram = TRUE)

methyl_proportion_cell <- function(orderObject, color = "yellow", 
                                makePlot=TRUE, ...){
  color <- tolower(color)
  if (color=="gch") color <- "yellow"
  color.indicator <- ifelse(color=="yellow", -1, 1)
  Proportion <- apply(orderObject$toClust, 1, function(x){
      sum(x == color.indicator * 3 | x == color.indicator * 4) / (length(x) / 2)
  })
  if (makePlot) {
    opar <- par(lwd=4)
    H = hist(Proportion, plot=FALSE, breaks = 15)
    plot(H, xlim=c(0,1), border=ifelse(color == "yellow", "gold2", "brown1"),
        col="gray75", main="Methylated Sites Per Cell/Molecule",
        lwd=2,...)
    par(opar)
  }
  return(Proportion)
}

#' Calculate the average methylation/accessibility status across all reads.
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
#' data(day7)
#' 
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")
#' methyl_average_status(orderObj, plotAverages = TRUE)
#' 
methyl_average_status <- function(orderObject, window_length = 20, 
                                makePlot = TRUE, ...)
{
    colLength <- ncol(orderObject$toClust)
    gch.num <- orderObject$toClust[,seq(1,(colLength / 2))]
    hcg.num <- orderObject$toClust[,seq((colLength / 2 + 1),colLength)]

    acc.sum <- colSums(gch.num == -3)
    meth.sum <- colSums(hcg.num == 3)
    acc.denom <- colSums(gch.num != 0)
    meth.denom <- colSums(hcg.num != 0)
    acc.avg <- acc.sum / acc.denom
    meth.avg <- meth.sum / meth.denom

    width <- window_length
    moving.acc.avg <- filter(x = acc.avg, filter = rep(1, width)) / width
    moving.meth.avg <- filter(x = meth.avg, filter = rep(1, width)) / width

    if (makePlot) {
        plot(moving.acc.avg, type = "l", col="gold2",
            xlab="Position along read", 
            ylab="Population-averaged status", ylim = c(0,1))
        lines(moving.meth.avg, col="brown1")
        legend("topright", legend=c("Methylation", "Accessibility"),
            fill=c("brown1", "gold2"))
    }
    return(list(meth_avg = moving.meth.avg, acc_avg = moving.acc.avg))
}
