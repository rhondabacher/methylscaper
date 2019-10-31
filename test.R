source("R/initialOrder.R")
source("R/refinement.R")
source("R/seqPlot.R")

load("data/day7.RData")
orderObj <- initialOrder(input.GCH = day7$gch, input.HCG = day7$hcg)

new.order <- refineFunction(orderObj, refineStart = 1, refineEnd = 25, Method = "PCA")

orderObj.new <- orderObj
orderObj.new$order1 <- new.order

plotSequence(orderObj)
plotSequence(orderObj.new)
