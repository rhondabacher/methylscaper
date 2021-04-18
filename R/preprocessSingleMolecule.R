#' Align the single-molecule data
#'
#' Runs the preprocessing methods for single-molecule data.
#'
#' @param ref A reference sequence to align the reads to.
#' @param fasta A list of reads/sequences from a single-molecule experiment (e.g. MAPit)
#' @param fasta_subset (optional) A vector of indices indicating which 
#'      sequences to process if a subset should be used. Leave this blank if
#'      all sequences should be processed.
#' @param multicoreParam (optional) A MulticoreParam object, used to align 
#'      sequences in parallel.
#' @param updateProgress (optional) Used to add a progress bar to the Shiny app. 
#'      Should not be used otherwise.
#' @param log_file (optional) String indicating where to save a log of the 
#'      alignment process. If left NULL, no log is saved. We highly recommend
#'      saving a log file.
#' @return The output is a list containing the the matrices 'gch' and 'hcg.
#'       Each is a dataframe with reads/cells on the rows and each column
#'      is a base-pair. The matrix represents the methylation state for cell
#'      across all base pairs. The coding is as follows:
#'          -2: unmethylated GCH or HCG site
#'          -1: base pairs between two unmethylated GCH or HCG sites
#'          0: base pairs between mismatching methylation states of 
#'              two GCH or HCG sites
#'          1: base pairs between two methylated GCH or HCG sites
#'          2: methylated GCH or HCG site
#' @importFrom Biostrings DNAString DNA_ALPHABET reverseComplement mismatchTable
#' @importFrom Biostrings pairwiseAlignment score alignedPattern alignedSubject
#' @importFrom seqinr c2s s2c read.fasta
#' @importFrom BiocParallel bplapply
#' @export
#' @examples 
#'  
#' data(reads_sm)
#' data(ref_seq)
#' example_alignedseq <- runAlign(fasta=reads_sm, ref = ref_seq[[1]], fasta_subset = 1:150)

runAlign <- function(ref, fasta, fasta_subset = seq(1,length(fasta)),
                     multicoreParam = NULL, updateProgress = NULL, 
                     log_file = NULL)
{
    fasta <- fasta[fasta_subset]
    ref_string <- DNAString(toupper(c2s(ref)))

    log_vector <- c("Beginning preprocessing")

    if (is.function(updateProgress)) {
      updateProgress(message = "Aligning sequences", value = 0.1)
    }
    alignment_out <- alignSequences(fasta, ref_string, log_vector,
                                      multicoreParam, updateProgress)

    alignedseq <- alignment_out$alignedseq
    log_vector <- alignment_out$log_vector

    if (is.function(updateProgress)) {
          updateProgress(message = "Identifying sites", value = 0.75)
    }
    # We want to avoid GCG sites:
    GCsites <- gregexpr("GC",c2s(ref_string),fixed=TRUE)[[1]] + 1
    CGsites <- gregexpr("CG",c2s(ref_string),fixed=TRUE)[[1]]

    cg_site_use <- s2c(paste(ref_string))[CGsites-1]
    gc_site_use <- s2c(paste(ref_string))[GCsites+1]

    log_vector <- c(log_vector, paste("Throwing out",
        length(which(gc_site_use == "G")) + length(which(cg_site_use == "G")), "GCG sites"))

    CGsites <- CGsites[which(cg_site_use != "G")]
    GCsites <- GCsites[which(gc_site_use != "G")]

    if (is.function(updateProgress)) {
      updateProgress(message = "Mapping sites", value = 0.8)
    }
     if (alignment_out$siteOrient == "rev"){
       GCsitesForMapping <- GCsites - 1
       CGsitesForMapping <- CGsites + 1
     } else {
       GCsitesForMapping <- GCsites
       CGsitesForMapping <- CGsites
     }

     if (is.null(multicoreParam)) {
       gcmap <- lapply(alignedseq, mapseq, 
                       sites=GCsitesForMapping, 
                       siteOrient=alignment_out$siteOrient)
       cgmap <- lapply(alignedseq, mapseq,
                       sites=CGsitesForMapping,
                       siteOrient=alignment_out$siteOrient)
     } else {
       gcmap <- bplapply(alignedseq, mapseq, 
                         sites=GCsitesForMapping,
                         siteOrient=alignment_out$siteOrient,
                         BPPARAM = multicoreParam)
       cgmap <- bplapply(alignedseq, mapseq, 
                         sites=CGsitesForMapping,
                         siteOrient=alignment_out$siteOrient,
                         BPPARAM = multicoreParam)
     }

    if (is.function(updateProgress)) {
          updateProgress(message = "Preparing matrices", value = 0.95)
    }

    saveCG <- data.matrix(do.call(rbind, lapply(cgmap, function(x) (x))))
    saveCG <- cbind(rownames(saveCG), saveCG)

    saveGC <- data.matrix(do.call(rbind, lapply(gcmap, function(x) (x))))
    saveGC <- cbind(rownames(saveGC), saveGC)
    if (is.function(updateProgress)) updateProgress(message = "Done", value = 1)

    return(list(hcg = saveCG, gch = saveGC, logs=log_vector))
}


# this handles the alignment of ALL the sequences, and returns the
# alignedseq object used in the runAlign function
# this needs the log_vector, multicoreParam, and updateProgress
# so that we can continue keeping track of these things
alignSequences <- function(fasta, ref_string, log_vector,
                        multicoreParam = NULL, updateProgress = NULL)
{
    ## this creates the substitution matrix for use in alignment
    penalty_mat <- matrix(0,length(DNA_ALPHABET[seq(1,4)]),
                            length(DNA_ALPHABET[seq(1,4)]))
    penalty_mat[seq(1,4),seq(1,4)] <- c(1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1)
    penalty_mat[penalty_mat==0] <- -5
    penalty_mat <- cbind(penalty_mat, c(0,0,0,0))
    penalty_mat <- rbind(penalty_mat, c(0,0,0,0,1))
    rownames(penalty_mat) <- colnames(penalty_mat) <- c(DNA_ALPHABET[seq(1,4)], "N")

    if (is.null(multicoreParam)) seqalign_out <- lapply(seq(1,length(fasta)),
        function(i) {
            if (is.function(updateProgress)) {
                updateProgress(message = "Aligning sequences",
                               detail = paste(i, "/", length(fasta)),
                               value = (0.1+ 0.65/length(fasta) * i))
            }
        seqalign(fasta[[i]], ref_string, substitutionMatrix = penalty_mat)
      }) else seqalign_out <- bplapply(seq(1,length(fasta)),
                                  function(i) seqalign(fasta[[i]], ref_string,
                                      substitutionMatrix = penalty_mat),
                                  BPPARAM = multicoreParam)
    useseqs <- lapply(seqalign_out, function(i) i$u)
    scores <- vapply(seqalign_out, function (i) i$score, numeric(1))
    maxAligns <- vapply(seqalign_out, function (i) i$maxAlign, numeric(1))

    score_cutoff_idx <- which.max(diff(sort(scores))) + 1
    score_cutoff <- sort(scores)[score_cutoff_idx]

    good_alignment_idxs <- which(scores > score_cutoff)

    if(length(good_alignment_idxs) == 0) {stop("No good alignments were found. See methylscaper FAQ for more details.")}

    alignedseq <- lapply(good_alignment_idxs, function(i){
        SEQ1 = s2c(paste(alignedPattern(useseqs[[i]])))
        SEQ2 = s2c(paste(alignedSubject(useseqs[[i]])))

        toreplace <- SEQ1[which(SEQ2=="-")]
        toreplace[toreplace!="C"] <- "."
        toreplace[toreplace=="C"] <- "."

        SEQ2[which(SEQ2=="-")] <- toreplace
        SEQ2 <- SEQ2[which(SEQ1!="-")]
        # if (maxAligns[i] == 1) SEQ2 <- s2c(paste(reverseComplement(DNAString(c2s(SEQ2)))))
        return(SEQ2)
    })
    log_vector <- c(log_vector, paste("Throwing out",
                    length(useseqs) - length(good_alignment_idxs), "alignments"))

    names(alignedseq) <- names(fasta)[good_alignment_idxs]
		
   getOrientation <- lapply(good_alignment_idxs, function(i){
     mismt = mismatchTable(useseqs[[i]])
     mismt = c(mismt$PatternSubstring, mismt$SubjectSubstring)
     mismm_nt <- names(sort(table(mismt), decreasing = T)[1:2])
     return(mismm_nt)
   })
   getOrientation <- table(unlist(getOrientation))
   if (!is.na(getOrientation["A"]) & !is.na(getOrientation["T"]) & getOrientation["A"] > getOrientation["T"]) {
     siteOrient <- "rev"
   } else if (is.na(getOrientation["C"]) & is.na(getOrientation["T"])) {
     siteOrient <- "rev"
   } else {siteOrient <- "stnd"}

   return(list(alignedseq = alignedseq, log_vector = log_vector, siteOrient = siteOrient))


}

# aligns a single read to the reference, returns the useseq string.
# Alignment is finished in the alignSequences fn
seqalign <- function(read, ref_string, substitutionMatrix) {

    fasta_string <- DNAString(toupper(c2s(read)))

    align_bb <- pairwiseAlignment(reverseComplement(ref_string),
                                 reverseComplement(fasta_string),
                                 type="global-local", gapOpening=8,
                                 substitutionMatrix=substitutionMatrix)
    align_ab <- pairwiseAlignment(ref_string, reverseComplement(fasta_string),
                                 type="global-local", gapOpening=8,
                                 substitutionMatrix=substitutionMatrix)
    align_aa <- pairwiseAlignment(ref_string, fasta_string,type="global-local",
                                 gapOpening=8,
                                 substitutionMatrix=substitutionMatrix)
    align_ba <- pairwiseAlignment(reverseComplement(ref_string), fasta_string,type="global-local",
                                 gapOpening=8,
                                 substitutionMatrix=substitutionMatrix)

    maxAlign <- which.max(c(score(align_aa), score(align_ab), score(align_bb), score(align_ba)))
    allseq <- list(align_aa, align_ab, align_bb, align_ba)
    useseq <- allseq[[maxAlign]]

    return(list(u = useseq, score = score(useseq), maxAlign = maxAlign))
}

mapseq <- function(i, sites, siteOrient) {
    editseq <- i

     if (siteOrient == "rev") {
       editseq[sites][editseq[sites] == "."] <- "A"
       editseq[sites][editseq[sites] == "A"] <- "-2"
       editseq[sites][editseq[sites] == "G"] <- "2"
       editseq[sites][editseq[sites] == "C"] <- "."
       editseq[sites][editseq[sites] == "T"] <- "."
     } else {
       editseq[sites][editseq[sites] == "."] <- "T"
       editseq[sites][editseq[sites] == "T"] <- "-2"
       editseq[sites][editseq[sites] == "C"] <- "2"
       editseq[sites][editseq[sites] == "G"] <- "."
       editseq[sites][editseq[sites] == "A"] <- "."
     }
     editseq[sites][editseq[sites] == "N"] <- "."
    # we need to make sure that the N sites stay marked with a "."
    missing_bp <- which(editseq == ".")

    sites_temp <- c(0, sites, length(editseq)+1)
    for (j in seq(1,(length(sites_temp)-1))) {
        tofill <- seq(sites_temp[j]+1,(sites_temp[j+1]-1))
        s1 <- editseq[pmax(1, sites_temp[j])]
        s2 <- editseq[pmin(length(i), sites_temp[j+1])]
        is_first <- (pmax(1, sites_temp[j]) == 1)

    if (s1 == "2" & s2 == "2") { fillvec <- 1 }
    else if (s1 == "2" & s2 == "-2") {fillvec <- 0}
    else if (s1 == "-2" & s2 == "2") {fillvec <- 0}
    else if (s1 == "-2" & s2 == "-2") {fillvec <- -1}
    else if (s1 == "." & s2 == "."){fillvec <- "."}
    else if (is_first & s2 == ".") {fillvec <- "."}
    else fillvec <- 0
    fillvec <- rep(fillvec, length(tofill))
    editseq[tofill] <- fillvec
    editseq[intersect(tofill, missing_bp)] <- "."
    }

    substring_table <- get_contig_substrings(editseq)
    long_missing <- which(substring_table$char == "." & substring_table$count > 3)
    short_non_missing <- which(substring_table$char == "-" &
                                      substring_table$count <= 20)
    short_non_missing_sr <- short_non_missing[(short_non_missing + 1)
                 %in% long_missing | (short_non_missing - 1) %in% long_missing]
    counts <- as.numeric(substring_table$count)
    # missing artifacts
    for (idx in short_non_missing_sr) {
      if (idx == 1) editseq[seq(1,counts[idx])] <- "."
      else {
          first <- sum(counts[seq(1,(idx-1))]) + 1
          last <- first + counts[idx] - 1
          editseq[first:last] <- "."
      }
    }
    return(editseq)
}


## we want to be able to get all contiguous substrings of a certain string...
# in particular one of the editseq strings used above
## i want to return a table
get_contig_substrings <- function(s) {
    ## we want some sort of table to keep track of the substrings
    substring_table <- data.frame(char = "a", count = 0)

    s <- ifelse(s == ".", ".", "-") ## use "-" to denote the non-missing ones,
    # so we only keep track of missing and non missing
    i <- 1
    while (i <= length(s)) {
        j <- i
        while(s[j] == s[i] & j <= length(s)) j <- j + 1
        substring_table <- rbind(substring_table, c(s[i], j - i))
        i <- j
    }
    substring_table <- data.frame(char = substring_table$char[-1],
                            count = as.numeric(substring_table$count[-1]))
    return(substring_table)
}
