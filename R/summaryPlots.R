# plot the percent of C at each methylation site
# this is either darkred or goldenyellow depending on the plot
percent_C <- function(){}

# the proportion of yellow at each read, histogram of proportions
proportion_yellow <- function(orderObject, plotHistogram=FALSE){
  proportion <- apply(orderObject$toClust, 1, function(x){
    sum(x == -3) / (length(x) / 2)
  })
  if (plotHistogram) hist(proportion, main="Proportion of Yellow per Read")
  return(proportion)
}