#' Process single-cell data
#'
#' This function subsets the data and prepares it for visualizing by
#' generating representation methylation-state matrices from
#' single-cell methylation data (for example, sc-MNT data).
#' We assume the data has already been preprocess using the subsetSC
#' function in methylscaper. See the vignette for a more thorough
#' explanation of each parameter. The output should be passed directly to the
#' plotting function.
#'
#' @param dataIn A list object containing two elements labelled gch and hcg (already pre-processed.)
#' @param startPos The index of the first position to include
#'  in the visualization. If using this within the R console it is
#'  recomended to specify the start and end directly.
#'  In the Shiny app, a slider will let the user refine these positions.
#' @param endPos The index of the final position to include in the visualization.
#' @param updateProgress A function for generating progress bars in the Shiny app.
#'   Should be left NULL otherwise.
#' @return The output is a list containing the elements 'gch' and 'hcg.
#'      Each is a dataframe with reads/cells on the rows and each column
#'      is a base-pair. The matrix is coded as follows:
#'          -2: unmethylated GCH or HCG site
#'          -1: base pairs between two unmethylated GCH or HCG sites
#'          0: base pairs between mismatching methylation states of
#'              two GCH or HCG sites
#'          1: base pairs between two methylated GCH or HCG sites
#'          2: methylated GCH or HCG site
#' @importFrom utils tail
#' @import data.table
#' @importFrom methods is
#' @export
#' @examples
#'
#' data(singlecell_subset)
#' prepsc.out <- prepSC(singlecell_subset,
#'     startPos = 105636488, endPos = 105636993
#' )
prepSC <- function(dataIn, startPos = NULL, endPos = NULL,
    updateProgress = NULL) {
    if (is(dataIn, "SummarizedExperiment") | is(dataIn, "SingleCellExperiment")) {
        dataIn <- reformatSCE(dataIn)
    }

    if (is.null(updateProgress) & (endPos - startPos) >= 50000) {
        message("Selected range is longer than 50k bp,
                    plot may take a few seconds to render.")
    }

    if (is.null(updateProgress) & (endPos - startPos) >= 100000) {
        message("Selected range is longer than 100k bp,
                    this may not be optimal.
                    Plot will render but may take a few minutes.")
    }

    gc_seq_data <- dataIn$gch
    cg_seq_data <- dataIn$hcg
    if (is.function(updateProgress)) {
        updateProgress(message = "Filtering HCG (endogenous methylation)
                                    site data", value = 0.1)
    }

    cg_seq_sub <- lapply(cg_seq_data, function(x) {
        tempSub <- x[order(x$pos), ]
        tempSub <- subset(tempSub, tempSub$pos >= startPos & tempSub$pos <= endPos)
        return(tempSub)
    })

    if (is.function(updateProgress)) {
        updateProgress(message = "Filtering GCH (accessibility/occupancy) site data", value = 0.5)
    }


    gc_seq_sub <- lapply(gc_seq_data, function(x) {
        tempSub <- x[order(x$pos), ]
        tempSub <- subset(tempSub, tempSub$pos >= startPos & tempSub$pos <= endPos)
        return(tempSub)
    })

    all_cg_sites <- unique(do.call(c, cg_seq_sub))
    all_gc_sites <- unique(do.call(c, gc_seq_sub))
    useseq <- intersect(
        which(vapply(cg_seq_sub, function(x) nrow(x), numeric(1)) > 1),
        which(vapply(gc_seq_sub, function(x) nrow(x), numeric(1)) > 1)
    )

    if (length(useseq) == 0) {
        return(NULL)
    }

    cg_seq_sub <- cg_seq_sub[useseq]
    gc_seq_sub <- gc_seq_sub[useseq]

    if (is.function(updateProgress)) {
        updateProgress(message = "Mapping HCG (endogenous methylation) site data", value = 0.75)
    }

    cg_outseq <- lapply(cg_seq_sub, function(x) {
        if (nrow(x) == 0) {
            NULL
        } else {
            mapSC(x, startPos, endPos)
        }
    })

    cg_outseq <- cg_outseq[!vapply(cg_outseq, is.null, logical(1))]
    hcg <- data.matrix(do.call(rbind, cg_outseq))
    rownames(hcg) <- as.character(seq(1, nrow(hcg)))

    if (is.function(updateProgress)) {
        updateProgress(message = "Mapping GCH (accessibility/occupancy) site data", value = 0.9)
    }

    gc_outseq <- lapply(gc_seq_sub, function(x) {
        if (nrow(x) == 0) {
            NULL
        } else {
            mapSC(x, startPos, endPos)
        }
    })

    gc_outseq <- gc_outseq[!vapply(gc_outseq, is.null, logical(1))]
    gch <- data.matrix(do.call(rbind, gc_outseq))
    rownames(gch) <- as.character(seq(1, nrow(gch)))

    return(list(gch = gch, hcg = hcg))
}

mapSC <- function(inSeq, startPos, endPos) {
    inSeq$basePos <- as.numeric(inSeq$pos) - startPos + 1
    fill_1 <- seq(startPos, endPos) - startPos + 1
    someMethyl <- which(inSeq$rate > 0)
    noMethyl <- which(inSeq$rate <= 0)
    fill_1[fill_1 %in% inSeq[someMethyl, ]$basePos] <- 2
    fill_1[fill_1 %in% inSeq[noMethyl, ]$basePos] <- -2
    fill_1[abs(fill_1) != 2] <- "."

    sites <- inSeq$basePos
    sites <- sites[sites > 0]
    editseq <- fill_1
    sites_temp <- c(0, sites, max(sites) + 1)

    for (j in seq(1, (length(sites_temp) - 1))) {
        if (sites_temp[j + 1] == 1) { # skip
        } else {
            tofill <- seq(sites_temp[j] + 1, (sites_temp[j + 1] - 1))
            s1 <- editseq[pmax(1, sites_temp[j])]
            s2 <- editseq[pmin(length(editseq), sites_temp[j + 1])]

            if (s1 == "2" & s2 == "2") {
                fillvec <- 1
            } else if (s1 == "2" & s2 == "-2") {
                fillvec <- 0
            } else if (s1 == "-2" & s2 == "2") {
                fillvec <- 0
            } else if (s1 == "-2" & s2 == "-2") {
                fillvec <- -1
            } else {
                fillvec <- 0
            }
            fillvec <- rep(fillvec, length(tofill))
            editseq[tofill] <- fillvec
        }
    }
    return(editseq)
}



#' Load in methylation data
#'
#' This function loads the single-cell files. It takes a path to the data files
#' and a chromosome number as arguments and returns the desired subset of the
#' data. Processing by chromosome is important for speed and memory efficiency.
#' The input files should be tab separated with three columns.
#' The first column is the chromosome, the second is the position (basepair), and the third
#' is the methylation indicator/rate. The folder should contain two subfolders titled
#' met and acc, with the endogenous methylation and accessibility methylation files,
#' respectively.
#'
#' @param path Path to the folder containing the single-cell files.
#' @param chromosome The chromosome to subset the files to.
#' @param startPos The index of the first position to include
#'  in the subsetting. This is optional as further narrowing of the
#'  position can be done in the visualization step/tab.
#'  In the Shiny app, a slider will let the user refine the positions.
#' @param endPos The index of the final position to include in subset.
#' @param updateProgress A function for generating progress bars in the Shiny app.
#'   Should be left NULL otherwise.
#' @return The output is RDS files that can be loaded into the visualization
#'  tab on the Shiny app
#' @import data.table
#' @export
#'
#' @examples
#' # example not run since needs directory input from user
#' # subsc.out <- subsetSC("filepath", chromosome=19)
subsetSC <- function(path, chromosome, startPos = NULL,
    endPos = NULL, updateProgress = NULL) {
    if (!is.list(path)) {
        cgfiles <- paste0(path, "/", "met/", sort(grep("met", list.files(paste0(path, "/met")), value = TRUE)))
        gcfiles <- paste0(path, "/", "acc/", sort(grep("acc", list.files(paste0(path, "/acc")), value = TRUE)))
    } else {
        cgfiles <- path[[1]]
        gcfiles <- path[[2]]

        getord <- order(path[[3]])
        cgfiles <- cgfiles[getord]

        getord <- order(path[[4]])
        gcfiles <- gcfiles[getord]
    }


    if (length(cgfiles) != length(gcfiles)) {
        stop("Must have the same number of methylation and
                accessibility files!")
    }

    useChr <- chromosome

    cg_seq <- list()
    for (i in seq(1, length(cgfiles))) {
        in_cg_seq <- fread(cgfiles[i],
            header = FALSE, stringsAsFactors = FALSE
        )
        # A couple possible formats from bismark:
        if (ncol(in_cg_seq) %in% c(4, 6)) {
            in_cg_seq <- in_cg_seq[, c(1, 2, 4)]
        }

        if (ncol(in_cg_seq) != 3) {
            stop("Data is not formatted correctly. See vignette for details.")
        }

        if (in_cg_seq[1, 1] == "chr") {
            in_cg_seq <- in_cg_seq[-1, ]
        }
        colnames(in_cg_seq) <- c("chr", "pos", "rate")
        in_cg_seq$pos <- as.numeric(in_cg_seq$pos)
        in_cg_seq$rate <- as.numeric(in_cg_seq$rate)

        in_cg_seq <- subset(in_cg_seq, in_cg_seq$chr == useChr)
        if (!is.null(startPos) & !is.null(endPos)) {
            in_cg_seq <- subset(in_cg_seq, in_cg_seq$pos >= startPos & in_cg_seq$pos <= endPos)
        }

        cg_seq[[i]] <- in_cg_seq
        if (is.function(updateProgress)) {
            updateProgress(message = "Reading HCG (endogenous methylation) site files", value = i / length(cgfiles))
        }
    }

    gc_seq <- list()
    for (i in seq(1, length(gcfiles))) {
        in_gc_seq <- fread(gcfiles[i], header = FALSE, stringsAsFactors = FALSE)
        if (ncol(in_gc_seq) %in% c(4, 6)) {
            in_gc_seq <- in_gc_seq[, c(1, 2, 4)]
        }
        if (ncol(in_gc_seq) != 3) {
            stop("Data is not formatted correctly. See vignette for details.")
        }

        colnames(in_gc_seq) <- c("chr", "pos", "rate")
        if (in_gc_seq[1, 1] == "chr") {
            in_gc_seq <- in_gc_seq[-1, ]
        }
        colnames(in_gc_seq) <- c("chr", "pos", "rate")
        in_gc_seq$pos <- as.numeric(in_gc_seq$pos)
        in_gc_seq$rate <- as.numeric(in_gc_seq$rate)

        in_gc_seq <- subset(in_gc_seq, in_gc_seq$chr == useChr)
        if (!is.null(startPos) & !is.null(endPos)) {
            in_gc_seq <- subset(in_gc_seq, in_gc_seq$pos >= startPos & in_gc_seq$pos <= endPos)
        }

        gc_seq[[i]] <- in_gc_seq
        if (is.function(updateProgress)) {
            updateProgress(message = "Reading GCH (accessibility/occupancy) site files", value = i / length(gcfiles))
        }
    }

    list(hcg = cg_seq, gch = gc_seq)
}
