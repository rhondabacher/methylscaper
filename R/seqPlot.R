#' Generate Sequence Plot
#'
#' Generates an ordered sequence plot of methylation data.
#'
#' @param orderObject An object of class \code{orderObject} that contains 
#'    the processed data and the ordering.
#' @param plotFast Logical, setting to FALSE will generate a higher 
#'    quality plot. TRUE generates a lower resolution file, useful to improve
#'    speed while testing. For publication quality use plotFast=TRUE.
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
#' data(singlemolecule_example)
#' 
#' orderObj <- initialOrder(singlemolecule_example, Method = "PCA")
#' plotSequence(orderObj)

plotSequence <- function(orderObject, plotFast=TRUE,
                           blankWidth=NULL, Title="",
                           drawLine=TRUE, drawKey=TRUE) {
    # Start with yellow at top as the default:
    toClust <- orderObject$toClust
    order1 <- orderObject$order1
    input_GCH <- toClust[,seq(1,(ncol(toClust) / 2))]
    input_HCG <- toClust[,seq((ncol(toClust) / 2 + 1),ncol(toClust))]

    myCols <- c("darkgoldenrod2", "yellow", "gray62", "black", 
                "gray80", "white", "gray80", "black", "gray62", "red", "darkred")
    myBreaks <- c(-5,-4,-3,-2.5,-2,-1,0,1,2,2.5,3,4)

    if (nrow(toClust) == 1) {
      input_HCG <- t(as.matrix(input_HCG))
      input_GCH <- t(as.matrix(input_GCH))
    }
    
    input_HCG_fix <- input_HCG
    input_GCH_fix <- input_GCH

		# This spacing often looks best:
		if (is.null(blankWidth)) blankWidth <- round(.12 * ncol(toClust) / 2) 
		
    blankCOLS <- matrix(rep(0, nrow(input_HCG)*blankWidth), 
                    nrow=nrow(input_HCG), ncol=blankWidth)
    toPlot_fix_og <- cbind(input_HCG_fix, blankCOLS, input_GCH_fix)

    sites = which(apply(abs(toPlot_fix_og), 2, function(x) any(x %in% c(4, 1))))
    sites[sites < (ncol(toPlot_fix_og)/2)] <- sites[sites < (ncol(toPlot_fix_og)/2)] - 1
    sites_scale <- sites / ncol(toPlot_fix_og) # relative to the 'plot'

    toPlot_fix <- toPlot_fix_og[rev(order1),]

    if (nrow(toClust) == 1) {
      toPlot_fix <- t(as.matrix(toPlot_fix))
    }
    ## Add some rows in order to add a legend:
    if(drawKey == TRUE) {
      totalHeight <- pmax(ceiling(.25*nrow(toPlot_fix) * 0.1), 3)
      blankROW <- matrix(rep(0, ncol(toPlot_fix)*totalHeight), 
                               nrow=totalHeight, ncol=ncol(toPlot_fix))
      keyHeight <- pmax(ceiling(totalHeight * 0.25), 2)
      blankROW[seq(1,keyHeight),seq(1,147)] <- 2 #147 is nucleosome
      blankROW[seq(1,keyHeight),(seq(1,147)+(ncol(input_HCG_fix) + ncol(blankCOLS)))] <- 2
      toPlot_fix <- rbind(blankROW,toPlot_fix)
    }
    # Plotting:
    par(xpd = FALSE, mar=c(2.1,2.2,4,1), mgp = c(0,.8,.1))
     image(t(toPlot_fix), col=myCols, axes=FALSE, breaks=myBreaks,
           useRaster=plotFast, ylim=c(0,1))
           
     # Where to put the axes:
     # convert these back to the site number so we can do refinement
     plot1 <- ncol(input_HCG_fix)/ncol(toPlot_fix) 
     plot2 <- (ncol(input_HCG_fix)+blankWidth)/ncol(toPlot_fix)
     
     sites = which(apply(toPlot_fix_og, 2, function(x) any(x %in% c(4, 1)))) - 1
     sites_scale <- sites / ncol(toPlot_fix_og) # relative to the 'plot'

     if (drawKey == TRUE) {
     axis(3, at = c(0,sites_scale, plot1), labels = FALSE,
          tick=TRUE, line = .9, col="black", cex=1, lwd=2, lwd.tick=0,
          col.ticks = "white", tck = -.02)
      }
     axis(3, at = c(sites_scale), labels = FALSE,
          tick=TRUE, line = .9, col=NA, cex=1,
          col.ticks = "black", tck = .02)

     sites = which(apply(toPlot_fix_og, 2, function(x) any(x %in% c(-4, -1))))
     sites_scale <- sites / ncol(toPlot_fix_og) # relative to the 'plot'

     if (drawKey == TRUE) {
     axis(3, at = c(plot2,sites_scale, 1), labels = FALSE,
          tick=TRUE, line = .9, col="black", cex=1, lwd=2, lwd.tick=0,
          col.ticks = "white", tck = -.02)
     }
     axis(3, at = c(sites_scale), labels = FALSE,
          tick=TRUE, line = .9, col=NA, cex=1,
          col.ticks = "black", tck = .02)
     title(Title, line=2.7)
  
     # these shifts of 0.025 seem arbitrary but tend to center the title a bit better
     title("HCG", adj = plot1 / 2 - 0.025, line = 1.5) 
     title("GCH", adj = plot2 + 0.025 + (1 - plot2) / 2, line = 1.5)
  
     toLabel <- rev(seq(1, length(order1), by=ceiling(length(order1)/8)))
     if (!(length(order1) %in% toLabel)) toLabel <- c(length(order1), toLabel)
     y_axis_starting_point <- ifelse(drawKey, nrow(blankROW) / nrow(toPlot_fix), 0)
     axis(2, at = seq(y_axis_starting_point,1,
                      length.out=length(toLabel)), labels=toLabel, lwd=2.5, cex.axis=1.3)
  
     toLabel <- round(c(seq(1, ncol(input_HCG_fix), length.out=5),
                        seq(1, ncol(input_GCH_fix), length.out=5)))
     axis(1, at = seq(0,plot1,length.out=5), labels=toLabel[seq(1,5)], lwd=2.5, cex.axis=1.3)
     axis(1, at = seq(plot2,1,length.out=5), labels=toLabel[seq(1,5)], lwd=2.5, cex.axis=1.3)
  }
