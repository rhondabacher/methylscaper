# plot the percent of C at each methylation site
# this is either darkred or goldenyellow depending on the plot
percent_C <- function(orderObject){
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
    plot(x = red.sites - ncol(dat)/2, y = c.red,
         col="red", pch=19, ylim=c(0,1), xlab = "Site Number", ylab="%C")
    lines(x = red.sites - ncol(dat)/2, y = c.red, col="red")
    points(x = yellow.sites, y = c.yellow, col="yellow", pch=19)
    lines(x = yellow.sites, y = c.yellow, col="yellow")
    n.sites <- length(union(red.sites, yellow.sites))
    labs <- union(red.sites, yellow.sites)[seq(1, n.sites, by=n.sites/12)]
    #axis(1, at=labs)
}

# the proportion of a given color at each read, histogram of proportions
# color defaults to yellow
proportion_color <- function(orderObject, color = "YELLOW", plotHistogram=FALSE){
  color.indicator <- ifelse(color=="YELLOW", -1, 1)
  Proportion <- apply(orderObject$toClust, 1, function(x){
    sum(x == color.indicator * 3 | x == color.indicator * 4) / (length(x) / 2)
  })
  if (plotHistogram) hist(Proportion,
                          main=paste("Proportion of", ifelse(color.indicator==-1, "Yellow", "Red"), "per Read"),
                          xlim=c(0,1))

  return(Proportion)
}

