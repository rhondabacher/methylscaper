seqalign <- function(i, fastq, ref.string) {
  
  fastq.string <- DNAString(toupper(c2s(fastq[[i]])))
  
  align.bb <- pairwiseAlignment(reverseComplement(ref.string), 
                                reverseComplement(fastq.string),type="global-local", gapOpening=8)
  align.ab <- pairwiseAlignment(ref.string, reverseComplement(fastq.string),type="global-local", gapOpening=8)
  align.aa <- pairwiseAlignment(ref.string, fastq.string,type="global-local", gapOpening=8)
  
  maxAlign <- which.max(c(score(align.bb), score(align.ab), score(align.aa)))
  allseq <- list(align.bb, align.ab, align.aa)
  useseq <- allseq[[maxAlign]]
  
  if (score(useseq) > -500) { # if all alignments are "bad" then throw them away.-500 worked well in practice.
    SEQ1 = s2c(paste(pattern(useseq)))
    SEQ2 = s2c(paste(subject(useseq)))
    toreplace <- SEQ1[which(SEQ2=="-")]
    toreplace[toreplace!="C"] <- "."
    toreplace[toreplace=="C"] <- "T"
    SEQ2[which(SEQ2=="-")] <- toreplace
    SEQ2 <- SEQ2[which(SEQ1!="-")]
    if (maxAlign == 1) {SEQ2 <- s2c(paste(reverseComplement(DNAString(c2s(SEQ2)))))}
    alignedseq <- SEQ2
  } else alignedseq <- NULL
  return(alignedseq)
}

mapseq <- function(i, sites) {
  editseq <- i
  editseq[sites][editseq[sites] == "T"] <- "-2"
  editseq[sites][editseq[sites] == "C"] <- "2"
  editseq[sites][editseq[sites] == "G"] <- "."
  editseq[sites][editseq[sites] == "A"] <- "."
  
  sites.temp <- c(0, sites, length(editseq)+1)  
  for (j in 1:(length(sites.temp)-1)) {
    tofill <- seq(sites.temp[j]+1,(sites.temp[j+1]-1))
    s1 <- editseq[pmax(1, sites.temp[j])]
    s2 <- editseq[pmin(650, sites.temp[j+1])]
    if (s1 == "2" & s2 == "2") {
      fillvec <- 1 } else if (s1 == "2" & s2 == "-2") {
        fillvec <- 0} else if (s1 == "-2" & s2 == "2") {
          fillvec <- 0} else if (s1 == "-2" & s2 == "-2") {
            fillvec <- -1} else {fillvec <- 0}
    fillvec <- rep(fillvec, length(tofill))
    editseq[tofill] <- fillvec
  }
  return(editseq)
}

#' @import Biostrings
#' @import seqinr
#' @import BiocParallel 
#' @export
runAlign <- function(ref, fastq, fastq.subset = (1:length(seq2)), multicoreParam = NULL, updateProgress = NULL)
{
  fastq <- fastq[fastq.subset]
  ref.string <- DNAString(toupper(c2s(ref[[1]])))
  
  penalty.mat <- matrix(0,length(DNA_ALPHABET[1:4]),length(DNA_ALPHABET[1:4]))
  penalty.mat[1:4,1:4] <- c(1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1)
  penalty.mat[penalty.mat==0] <- -2
  rownames(penalty.mat) <- colnames(penalty.mat) <- DNA_ALPHABET[1:4]
  
  if (is.function(updateProgress)) updateProgress(message = "Aligning sequences", value = 0.1)
  
  if (is.null(multicoreParam)) alignedseq <- lapply(1:length(fastq), function(i) {
    if (is.function(updateProgress))updateProgress(message = "Aligning seqences",
                                                                 detail = paste(i, "/", length(fastq)), 
                                                                 value = (0.1+ 0.65/length(fastq) * i))
    seqalign(i, fastq, ref.string)
    })
  else alignedseq <- bplapply(1:length(fastq), function(i) seqalign(i, fastq, ref.string), BPPARAM = multicoreParam)
  names(alignedseq) <- names(fastq)
  
  # Only keep the 'good' alignments
  alignedseq <- alignedseq[which(!sapply(alignedseq, is.null))]
  
  if (is.function(updateProgress)) updateProgress(message = "Identifying sites", value = 0.75)
  # We want to avoid GCG sites:
  GCsites <- gregexpr("GC",c2s(ref.string),fixed=TRUE)[[1]] + 1
  GCsites <- GCsites[which(s2c(paste(ref.string))[GCsites+1] != "G")]
  
  CGsites <- gregexpr("CG",c2s(ref.string),fixed=TRUE)[[1]]
  CGsites <- CGsites[which(s2c(paste(ref.string))[CGsites-1] != "G")]
  
  if (is.function(updateProgress)) updateProgress(message = "Mapping sites", value = 0.8)
  
  if (is.null(multicoreParam))
  {
    gcmap <- lapply(alignedseq, mapseq, sites=GCsites)
    cgmap <- lapply(alignedseq, mapseq, sites=CGsites)
  }
  else
  {
    gcmap <- bplapply(alignedseq, mapseq, sites=GCsites, BPPARAM = multicoreParam)
    cgmap <- bplapply(alignedseq, mapseq, sites=CGsites, BPPARAM = multicoreParam)
  }
  
  if (is.function(updateProgress)) updateProgress(message = "Preparing matrices", value = 0.95)
  
  saveCG <- data.matrix(do.call(rbind, lapply(cgmap, function(x) (x))))
  saveCG <- cbind(rownames(saveCG), saveCG)
  
  saveGC <- data.matrix(do.call(rbind, lapply(gcmap, function(x) (x))))
  saveGC <- cbind(rownames(saveGC), saveGC)
  if (is.function(updateProgress)) updateProgress(message = "Done", value = 1)
  return(list(hcg = saveCG, gch = saveGC))
}
