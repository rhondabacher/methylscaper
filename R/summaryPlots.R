#' Calculates the percentage of C at each methylation site
#'
#' @param orderObject An object of class \code{orderObject}
#' @param plotPercents Logical, indicates whether to generate the percentage plot
#' @param ... Additional parameters used by the \code{plot} function.
#'
#' @importFrom graphics hist lines plot points
#' @export
percent_C <- function(orderObject, plotPercents = FALSE, ...){
    dat <- orderObject$toClust
    red.sites <- which(dat[1,1:ncol(dat)] == 4 |
                          dat[1,1:ncol(dat)] == 1)

    yellow.sites <- which(dat[1,1:ncol(dat)] == -4 |
                          dat[1,1:ncol(dat)] == -1)
    c.red <- sapply(red.sites, function(i){
        sum(dat[, i] == 4) / nrow(dat)
    })
    c.yellow <- sapply(yellow.sites, function(i){
        sum(dat[, i] == -4) / nrow(dat)
    })
    if (plotPercents)
    {
       plot(x = red.sites - ncol(dat)/2, y = c.red,
            col="brown1", pch=19, ylim=c(0,1), xlab = "Region (Base Pair)", ylab="%C",
            bty='n', cex.lab=1.3, xaxt='n', yaxt='n',...)
       axis(side = 1, lwd = 2, cex.axis=1.2)
       axis(side = 2, lwd = 2, cex.axis=1.2)
       lines(x = red.sites - ncol(dat)/2, y = c.red, col="brown1")
       points(x = yellow.sites, y = c.yellow, col="gold2", pch=19)
       lines(x = yellow.sites, y = c.yellow, col="gold2")

       n.sites <- length(union(red.sites, yellow.sites))
       labs <- union(red.sites, yellow.sites)[seq(1, n.sites, by=n.sites/12)]
    }
    final <- list(c.red, c.yellow)
    names(final) = c("red", "yellow")
    return(final)

}

#' Calculate the proportion of methylated bases for the GCH and HCG data sets.
#'
#' @param orderObject An object of class \code{orderObject}
#' @param color Indicates which data set to compute proportions for
#' @param plotHistogram Indicates whether to plot a histogram of the proportions across all reads.
#' @param ... Additional parameters used by the \code{hist} function.
#'
#' @importFrom graphics hist
#' @export
proportion_color <- function(orderObject, color = "YELLOW", plotHistogram=FALSE, ...){
  color.indicator <- ifelse(color=="YELLOW", -1, 1)
  Proportion <- apply(orderObject$toClust, 1, function(x){
    sum(x == color.indicator * 3 | x == color.indicator * 4) / (length(x) / 2)
  })
  if (plotHistogram) {
    opar <- par(lwd=4)
    H = hist(Proportion, plot=F)
    plot(H, xlim=c(0,1), border=ifelse(color == "YELLOW", "gold2", "brown1"),
       col="gray75",
       lwd=2,...)
    par(opar)
  }
  return(Proportion)
}

#' Calculate the average methylation/accessibility status across all reads.
#'
#' @param orderObject An object of class \code{orderObject}
#' @param window_length Length of the window to be used to compute a moving average.
#' @param plotAverages Logical, indicates whether to generate a line plot of average status.
#' @param ... Addition parameters used by the \code{plot} function.
#'
#' @importFrom stats filter
#' @importFrom graphics legend
#' @export
average_status <- function(orderObject, window_length = 1, plotAverages = FALSE, ...)
{
    gch.num <- orderObject$toClust[,1:(ncol(orderObject$toClust) / 2)]
    hcg.num <- orderObject$toClust[,(ncol(orderObject$toClust) / 2 + 1):ncol(orderObject$toClust)]

    acc.sum <- colSums(gch.num == -3)
    meth.sum <- colSums(hcg.num == 3)
    acc.denom <- colSums(gch.num != 0)
    meth.denom <- colSums(hcg.num != 0)
    acc.avg <- acc.sum / acc.denom
    meth.avg <- meth.sum / meth.denom

    width <- window_length
    moving.acc.avg <- filter(x = acc.avg, filter = rep(1, width)) / width
    moving.meth.avg <- filter(x = meth.avg, filter = rep(1, width)) / width

    if (plotAverages)
    {
        plot(moving.acc.avg, type = "l", col="gold2",
             xlab="Position along read", ylab="Population-averaged status", ylim = c(0,1))
        lines(moving.meth.avg, col="brown1")
        legend("topright", legend=c("Methylation", "Accessibility"), fill=c("brown1", "gold2"))
    }
    return(list(meth_avg = moving.meth.avg, acc_avg = moving.acc.avg))
}
