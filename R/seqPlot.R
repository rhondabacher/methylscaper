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
#' @param shinySizer internal sizing parameter for plot layout in shiny app.
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
plotSequence <- function(
        orderObject, plotFast = TRUE,
        blankWidth = NULL, Title = "",
        drawLine = TRUE, drawKey = TRUE, shinySizer = 0) {
    # Start with yellow at top as the default:
    toClust <- orderObject$toClust
    order1 <- orderObject$order1
    input_GCH <- toClust[, seq(1, (ncol(toClust) / 2))]
    input_HCG <- toClust[, seq((ncol(toClust) / 2 + 1), ncol(toClust))]

    myCols <- c(
        "darkgoldenrod2", "yellow", "gray62", "black",
        "gray80", "white", "gray80", "black", "gray62", "red", "darkred"
    )
    myBreaks <- c(-5, -4, -3, -2.5, -2, -1, 0, 1, 2, 2.5, 3, 4)

    if (nrow(toClust) == 1) {
        input_HCG <- t(as.matrix(input_HCG))
        input_GCH <- t(as.matrix(input_GCH))
    }

    input_HCG_fix <- input_HCG
    input_GCH_fix <- input_GCH

    # This spacing often looks best:
    if (is.null(blankWidth)) blankWidth <- round(.12 * ncol(toClust) / 2)

    blankCOLS <- matrix(rep(0, nrow(input_HCG) * blankWidth),
        nrow = nrow(input_HCG), ncol = blankWidth
    )
    toPlot_fix_og <- cbind(input_HCG_fix, blankCOLS, input_GCH_fix)

    sites <- which(apply(abs(toPlot_fix_og), 2, function(x) any(x %in% c(4, 1))))
    sites_scale <- sites / ncol(toPlot_fix_og) # relative to the 'plot'

    toPlot_fix <- toPlot_fix_og[rev(order1), ]

    if (nrow(toClust) == 1) {
        toPlot_fix <- t(as.matrix(toPlot_fix))
    }
    ## Add some rows in order to add a legend:
    if (drawKey == TRUE) {
        totalHeight <- pmax(ceiling(nrow(toPlot_fix) * 0.1), 3)
        blankROW <- matrix(rep(0, ncol(toPlot_fix) * totalHeight),
            nrow = totalHeight, ncol = ncol(toPlot_fix)
        )
        keyHeight <- pmax(ceiling(totalHeight * 0.5), 2)
        blankROW[seq(1, keyHeight), seq(1, 147)] <- 2 # 147 is nucleosome
        blankROW[seq(1, keyHeight), (seq(1, 147) + (ncol(input_HCG_fix) + ncol(blankCOLS)))] <- 2
        toPlot_fix <- rbind(blankROW, toPlot_fix)
    }
    # Plotting:
    if (shinySizer >= 750) TCK <- -.015 else TCK <- -.02
    if (shinySizer >= 750) TLINE <- .6 else TLINE <- .5
    if (shinySizer >= 750) topMAR <- 6 else topMAR <- 4

    par(xpd = FALSE, mar = c(2, 2, topMAR, 1), mgp = c(0, 0.5, 0))
    image(t(toPlot_fix),
        col = myCols, axes = FALSE, breaks = myBreaks,
        useRaster = plotFast, ylim = c(0, 1)
    )
    axis(3,
        at = sites_scale, labels = rep("", length(sites_scale)),
        tick = TRUE, line = TLINE, col = "white", cex = 1, lwd = 1,
        col.ticks = "black", tck = TCK
    )
    title(Title, line = 2, col.main = "royalblue4")

    # Where to put the axes:
    # convert these back to the site number so we can do refinement
    plot1 <- round(ncol(input_HCG_fix) / ncol(toPlot_fix), 2)
    plot2 <- round((ncol(input_HCG_fix) + blankWidth) / ncol(toPlot_fix), 2)

    if (shinySizer >= 750) lSHIFT <- 2.2 else lSHIFT <- 1.5
    # these shifts of 0.025 seem arbitrary but tend to center the title a bit better
    title("HCG", adj = plot1 / 2 - 0.025, line = lSHIFT)
    title("GCH", adj = plot2 + 0.025 + (1 - plot2) / 2, line = lSHIFT)

    toLabel <- rev(seq(1, length(order1), by = ceiling(length(order1) / 8)))
    if (!(length(order1) %in% toLabel)) toLabel <- c(length(order1), toLabel)
    y_axis_starting_point <- ifelse(drawKey, nrow(blankROW) / nrow(toPlot_fix), 0)
    axis(2, at = seq(y_axis_starting_point, 1,
        length.out = length(toLabel)
    ), labels = toLabel)

    toLabel <- round(c(
        seq(1, ncol(input_HCG_fix), length.out = 5),
        seq(1, ncol(input_GCH_fix), length.out = 5)
    ))
    axis(1, at = seq(0, plot1, length.out = 5), labels = toLabel[seq(1, 5)])
    axis(1, at = seq(plot2, 1, length.out = 5), labels = toLabel[seq(1, 5)])

    # Just drawing a straight DNA line:
    if (drawLine == TRUE) {
        par(xpd = NA)
        top1 <- 1.01
        segments(0, top1, plot1, top1, lwd = 1)
        segments(plot2, top1, 1, top1, lwd = 1)
    }
}
