seqalign <- function(i) {
  seq2string <- DNAString(toupper(c2s(seq2[[i]])))
  
  align.bb <- pairwiseAlignment(reverseComplement(seq1string), reverseComplement(seq2string),type="global-local", gapOpening=8)
  align.ab <- pairwiseAlignment(seq1string, reverseComplement(seq2string),type="global-local", gapOpening=8)
  align.aa <- pairwiseAlignment(seq1string, seq2string,type="global-local", gapOpening=8)
  
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


runAlign <- function(seq1.file, seq2.file, gch.file.name, hcg.file.name)
{
  seq1 <- read.fasta(seq1.file)
  seq2 <- read.fasta(seq2.file)
  seq1string <- DNAString(toupper(c2s(seq1[[1]])))
  
  penalty.mat <- matrix(0,length(DNA_ALPHABET[1:4]),length(DNA_ALPHABET[1:4]))
  penalty.mat[1:4,1:4] <- c(1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1)
  penalty.mat[penalty.mat==0] <- -2
  rownames(penalty.mat) <- colnames(penalty.mat) <- DNA_ALPHABET[1:4]
  
  alignedseq <- bplapply(1:length(seq2), seqalign) # this is an apply function to handle parallel processing
  names(alignedseq) <- names(seq2)
  
  # Only keep the 'good' alignments
  alignedseq <- alignedseq[which(!sapply(alignedseq, is.null))]
  
  # We want to avoid GCG sites:
  GCsites <- gregexpr("GC",c2s(seq1string),fixed=TRUE)[[1]] + 1
  GCsites <- GCsites[which(s2c(paste(seq1string))[GCsites+1] != "G")]
  
  CGsites <- gregexpr("CG",c2s(seq1string),fixed=TRUE)[[1]]
  CGsites <- CGsites[which(s2c(paste(seq1string))[CGsites-1] != "G")]
  
  gcmap <- bplapply(alignedseq, mapseq, sites=GCsites)
  cgmap <- bplapply(alignedseq, mapseq, sites=CGsites)
  
  saveCG <- data.matrix(do.call(rbind, lapply(cgmap, function(x) (x))))
  saveCG <- cbind(rownames(saveCG), saveCG)
  write.table(saveCG, file=hcg.file.name, quote=F, row.names = F, sep="\t")
  
  saveGC <- data.matrix(do.call(rbind, lapply(gcmap, function(x) (x))))
  saveGC <- cbind(rownames(saveGC), saveGC)
  write.table(saveGC, file=gch.file.name, quote=F, row.names = F, sep="\t")
}
