#' Generate Sequence Plot
#'
#' Generates an ordered sequence plot of methylation data.
#'
#' @param orderObject An object of class \code{orderObject} that contains 
#'    the processed data and the ordering.
#' @param plotFAST Logical, setting to FALSE will generate a higher 
#'    quality plot. TRUE generates a lower resolution file, useful to improve
#'    speed while testing. For publication quality use plotFAST=TRUE.
#' @param blankWidth Indicates the amount of space to leave between
#'     the two plots
#' @param Title The title of the plot.
#' @param drawLine Logical, indicates whether to draw a line above 
#'    the CG/GC sites.
#' @param drawKey Logical, indicates whether to draw a key representing
#'    a 147bp nucleosome at the bottom of the plot.
#'
#' @return Output is two side-by-side heatmaps with the endogenous 
#'    methylation (HCG) on the left and the acciessibility methylation (GCH)
#'    on the right. The tick marks at the top indicate
#'    either HCG or GCH sites when drawLine=TRUE. If drawKey=TRUE then a black
#'    rectangle key is plot at the bottom of the heatmap that is 147 basepairs
#'    long. In the HCG plot, red patches represent methylation between two 
#'    sites; black patches represent unmethylated bases between two
#'    unmethylated sites; and gray patches are base pairs which have one
#'    methylated site and one unmethylated site flanking. In the HCG plot,
#'    yellow patches represent accessibility between two sites; black patches
#'    represent occupied bases  between two occupied sites; and gray patches
#'    are base pairs which have one methylated site and one 
#'    unmethylated site flanking. 
#'
#' @importFrom graphics abline image par axis segments title
#' @export
#' @examples 
#'  
#' data(day7)
#' 
#' orderObj <- initialOrder(day7$gch, day7$hcg, Method = "PCA")
#' plotSequence(orderObj)

plotSequence <- function(orderObject, plotFAST=TRUE,
                           blankWidth=75, Title="",
                           drawLine=TRUE, drawKey=TRUE) {
    # Start with yellow at top as the default:
    toClust <- orderObject$toClust
    order1 <- orderObject$order1
    input_GCH <- toClust[,seq(1,(ncol(toClust) / 2))]
    input_HCG <- toClust[,seq((ncol(toClust) / 2 + 1),ncol(toClust))]

    mycols <- c("darkgoldenrod2", "yellow", "gray62", "black", 
                "gray80", "white", "gray80", "black", "gray62", "red", "darkred")
    VALS <- c(-5,-4,-3,-2.5,-2,-1,0,1,2,2.5,3,4)

    input_HCG_fix <- input_HCG
    input_GCH_fix <- input_GCH

    blankCOLS <- matrix(rep(0, nrow(input_HCG)*blankWidth), 
                    nrow=nrow(input_HCG), ncol=blankWidth)
    toPlot_fix_og <- cbind(input_HCG_fix, blankCOLS, input_GCH_fix)

    sites = which(apply(abs(toPlot_fix_og), 2, function(x) any(x %in% c(4, 1))))
    sites_scale <- sites / ncol(toPlot_fix_og) # relative to the 'plot'

    toPlot_fix <- toPlot_fix_og[rev(order1),]

    ## Add some rows in order to add a legend:
    if(drawKey == TRUE) {
      blankROW <- matrix(rep(0, ncol(toPlot_fix)*12), 
                          nrow=12, ncol=ncol(toPlot_fix))
      blankROW[seq(5,8),seq(1,147)] <- 2 #147 is nucleosome
      blankROW[seq(5,8),seq((seq(1,147)+(ncol(input_HCG_fix) + ncol(blankCOLS))))] <- 2
      toPlot_fix <- rbind(blankROW,toPlot_fix)
    }
    # Plotting:
    par(xpd = FALSE, mar=c(2,2,2,1), mgp = c(0,0.5,0))
    image(t(toPlot_fix), col=mycols, axes=FALSE, breaks=VALS,
        main=Title, useRaster=plotFAST, ylim=c(0,1.028))
    axis(3, at = sites_scale, labels = rep("", length(sites_scale)),
        tick=TRUE, line = .1, col="white", cex=1, lwd=1,
        col.ticks = "black", tck = .02)
    # Where to put the axes:
    # convert these back to the site number so we can do refinement
    plot1 <- round(ncol(input_HCG_fix)/ncol(toPlot_fix), 2) 
    plot2 <- round((ncol(input_HCG_fix)+blankWidth)/ncol(toPlot_fix), 2)

    # these shifts of 0.025 seem arbitrary but tend to center the title a bit better
    title("HCG", adj = plot1 / 2 - 0.025) 
    title("GCH", adj = plot2 + 0.025 + (1 - plot2) / 2)

    toLabel <- rev(seq(1, length(order1), by=round(length(order1)/8)))
    if (!(length(order1) %in% toLabel)) toLabel <- c(length(order1), toLabel)
    y_axis_starting_point <- ifelse(drawKey, nrow(blankROW) / nrow(toPlot_fix), 0)
    axis(2, at = seq(y_axis_starting_point,1,
                    length.out=length(toLabel)), labels=toLabel)

    toLabel <- round(c(seq(1, ncol(input_HCG_fix), length.out=5),
                        seq(1, ncol(input_GCH_fix), length.out=5)))
    axis(1, at = seq(0,plot1,length.out=5), labels=toLabel[seq(1,5)])
    axis(1, at = seq(plot2,1,length.out=5), labels=toLabel[seq(1,5)])

    # Just drawing a straight DNA line:
    if (drawLine==TRUE) {
        par(xpd=NA)
        top1 <- 1.01
        segments(0,top1,plot1,top1, lwd=1)
        segments(plot2,top1, 1,top1, lwd=1)
    }
  }
