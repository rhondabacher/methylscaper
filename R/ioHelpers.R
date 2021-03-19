#' This is an internal function for now and a place-holder 
#' in case SingleCellExperiment may be
#' used in the future. We may need to update this later. Assumes 
#' rownames are formatted like chr_pos and there are two assays.
#' The assay names are `methylation_met` for endogenous methylation 
#' and `methylation_acc` for accessibility.
#' @param dataIn Input SCE object passed to the prepSC() function.
#' @return List object containing the met and acc tables.
#' @importFrom SummarizedExperiment assay
reformatSCE <- function(dataIn) {
  
   check_package("SingleCellExperiment")
	 check_package("SummarizedExperiment")
  
   dataOut.gch <- assay(dataIn, "methylation_acc")
   locs_list <- rownames(dataOut.gch)
   locs_list <- data.frame(do.call(rbind, strsplit(locs_list, "_")))
   colnames(locs_list) <- c("chr", "pos")
   locs_list$pos <- as.numeric(locs_list$pos)
   dataOut.gch.list <- split(dataOut.gch, rep(seq(1,ncol(dataOut.gch)), each = nrow(dataOut.gch)))
	 
   
   dataOut.hcg <- assay(dataIn, "methylation_met")
   locs_list <- rownames(dataOut.hcg)
   locs_list <- data.frame(do.call(rbind, strsplit(locs_list, "_")))
   colnames(locs_list) <- c("chr", "pos")
   locs_list$pos <- as.numeric(locs_list$pos)
   dataOut.hcg.list <- split(dataOut.hcg, rep(seq(1,ncol(dataOut.hcg)), each = nrow(dataOut.hcg)))
   # dataOut.hcg.list <- lapply(dataOut.hcg.list, function(x) GRanges(seqnames = locs_list$chr,
   #                                         ranges = IRanges(start = locs_list$pos, width = 1L),
   #                                         rate=x))

   return(list(gch=dataOut.gch.list, hcg=dataOut.hcg.list))
  }


#' Check if a package is installed
#' @param package Name of the package.
#' @return Message is returned if package not installed.
check_package <- function(package){
  suppressMessages({
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(package," package is needed, please install it.")
    }
  })
}
