#' Calculates the percentage of C at each methylation site
#'
#' @param orderObject An object of class \code{orderObject}
#' @param plotPercents Logical, indicates whether to generate the percentage plot
#' @param ... Additional parameters used by the \code{plot} function.
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
             col="darkred", pch=19, ylim=c(0,1), xlab = "Region (Base Pair)", ylab="%C", ...)
        lines(x = red.sites - ncol(dat)/2, y = c.red, col="darkred")
        points(x = yellow.sites, y = c.yellow, col="goldenrod", pch=19)
        lines(x = yellow.sites, y = c.yellow, col="goldenrod")
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
#' @export
proportion_color <- function(orderObject, color = "YELLOW", plotHistogram=FALSE, ...){
  color.indicator <- ifelse(color=="YELLOW", -1, 1)
  Proportion <- apply(orderObject$toClust, 1, function(x){
    sum(x == color.indicator * 3 | x == color.indicator * 4) / (length(x) / 2)
  })
  if (plotHistogram) hist(Proportion,
                          xlim=c(0,1),
                          border=ifelse(color == "YELLOW", "goldenrod", "darkred"),
                          col="gray",
                          lwd=2, ...)

  return(Proportion)
}

