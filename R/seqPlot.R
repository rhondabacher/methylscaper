#Plotting customization:
  ## blankWidth: How much space between the two plots
  ## drawLine: Whether a line is drawn for the "DNA" strand is drawn above the CG/GC sites.
  ## Title: Title of plot
  ## plotFAST: fast renders a low quality version. for high quality publication image set FALSE.
  ## drawKey: whether to draw a key representing a nucleosome of 147bp at the bottom.

#' Generate Sequence Plot
#' 
#' Generates an ordered sequence plot of methylation data.
#' 
#' @param orderObject An object of class \code{orderObject} that contains the processed data and the ordering.
#' @param plotFAST Logical, setting to FALSE will generate a higher quality plot.
#' @param blankWidth Indicates the amount of space to leave between the two plots
#' @param Title The title of the plot.
#' @param drawLine Logical, indicates whether to draw a line above the CG/GC sites.
#' @param drawKey Logical, indicates whether to draw a key representing a 147bp nucleosome at the bottom of the plot.
#' 
#' @importFrom graphics abline image par axis segments
#' @export
plotSequence <- function(orderObject, plotFAST=TRUE,
                           blankWidth=150, Title="",
                           drawLine=T, drawKey=T) {
    # Start with yellow at top as the default:
    toClust <- orderObject$toClust
    order1 <- orderObject$order1
    input.GCH <- toClust[,1:(ncol(toClust) / 2)]
    input.HCG <- toClust[,(ncol(toClust) / 2 + 1):ncol(toClust)]
    
    mycols <- c("darkgoldenrod2", "yellow", "gray62",  "black", "gray80", "white", "gray80", "black", "gray62", "red", "darkred")
    VALS <- c(-5,-4,-3,-2.5,-2,-1,0,1,2,2.5,3,4)
    
    
    input.HCG.fix <- input.HCG
    input.GCH.fix <- input.GCH
    
    blankCOLS <- matrix(rep(0, nrow(input.HCG)*blankWidth), nrow=nrow(input.HCG), ncol=blankWidth)
    toPlot.fix.og <- cbind(input.HCG.fix, blankCOLS, input.GCH.fix)
    
    sites = which(apply(abs(toPlot.fix.og), 2, function(x) any(x %in% c(4, 1))))
    sites.scale <- sites / ncol(toPlot.fix.og) # relative to the 'plot'
    
    
    toPlot.fix <- toPlot.fix.og[rev(order1),]
    
    ## Add some rows in order to add a legend:
    if(drawKey == TRUE) {
      blankROW <- matrix(rep(0, ncol(toPlot.fix)*12), nrow=12, ncol=ncol(toPlot.fix))
      blankROW[5:8,1:147] <- 2
      blankROW[5:8,((1:147)+(ncol(input.HCG.fix) + ncol(blankCOLS)))] <- 2
      toPlot.fix <- rbind(blankROW,toPlot.fix)
    }
    
    # Plotting:
    par(xpd = F, mar=c(2,2,2,1))
    image(t(toPlot.fix), col=mycols, axes=F, breaks=VALS, 
          main=Title, useRaster=plotFAST, ylim=c(0,1.03))
    axis(3, at = sites.scale, labels = rep("", length(sites.scale)),
         tick=T, line = .1, col="white", cex=1, lwd=1.5,
         col.ticks = "black", tck = .02)
    # Where to put the axes:
    plot1 <- round(ncol(input.HCG.fix)/ncol(toPlot.fix), 2) # convert these back to the site number so we can do refinement
    plot2 <- round((ncol(input.HCG.fix)+blankWidth)/ncol(toPlot.fix), 2)
    
    toLabel <- rev(c(seq(1, length(order1), by=round(length(order1)/8)), length(order1)))
    axis(2, at = seq(0.077,1,length.out=length(toLabel)), labels=toLabel)
    toLabel <- round(c(seq(1, ncol(input.HCG.fix), length.out=5),
                       seq(1, ncol(input.GCH.fix), length.out=5)))
    axis(1, at = c(seq(0,plot1,length.out=5), seq(plot2,1,length.out=5)), labels=toLabel)
    
    # Just drawing a straight DNA line:
    if (drawLine==TRUE) {
      par(xpd=NA)
      top1 <- 1.01
      segments(0,top1,plot1,top1, lwd=1)
      segments(plot2,top1, 1,top1, lwd=1)
    } 
  }
