#' Process single-cell data
#'
#' @export
prepSC <- function(gc_seq_data, cg_seq_data, startPos, endPos,
                   updateProgress = NULL)
{
    if (is.function(updateProgress))
        updateProgress(message = "Filtering CG data", value = 0.1)
  cg_seq_sub <- lapply(cg_seq_data, function(x) {
    QQ <- x[order(x$pos),]
    QQ = subset(QQ, pos >= startPos & pos <= endPos)
    return(QQ)
  })

    if (is.function(updateProgress))
        updateProgress(message = "Filtering GC data", value = 0.5)
  gc_seq_sub <- lapply(gc_seq_data, function(x) {
    QQ <- x[order(x$pos),]
    QQ = subset(QQ, pos >= startPos & pos <= endPos)
    return(QQ)
  })

  all_cg_sites <- unique(do.call(c, cg_seq_sub))
  all_gc_sites <- unique(do.call(c, gc_seq_sub))

  useseq <- intersect(which(sapply(cg_seq_sub, function(x) nrow(x)) > 0 ),
                      which(sapply(gc_seq_sub, function(x) nrow(x)) > 0 ))
  cg_seq_sub <- cg_seq_sub[useseq]
    gc_seq_sub <- gc_seq_sub[useseq]


    if (is.function(updateProgress))
        updateProgress(message = "Mapping CG data", value = 0.75)
  cg_outseq <- lapply(cg_seq_sub, function(x) mapSC(x, startPos, endPos))
  hcg <- data.matrix(do.call(rbind, cg_outseq))
  rownames(hcg) <- as.character(1:nrow(hcg))

    if (is.function(updateProgress))
        updateProgress(message = "Mapping GC data", value = 0.9)
  gc_outseq <- lapply(gc_seq_sub, function(x) mapSC(x, startPos, endPos))
  gch <- data.matrix(do.call(rbind, gc_outseq))
  rownames(gch) <- as.character(1:nrow(gch))

  list(gch = gch, hcg = hcg)

}

mapSC <- function(IN.seq, startPos, endPos) {
  IN.seq$pos <- IN.seq$pos - startPos + 1
  fill.1 <- seq(startPos, endPos) - startPos + 1
  someMethyl <- which(IN.seq$rate > 0)
  noMethyl <- which(IN.seq$rate <= 0)
  fill.1[fill.1 %in% IN.seq[someMethyl,]$pos] <- 2
  fill.1[fill.1 %in% IN.seq[noMethyl,]$pos] <- -2
  fill.1[abs(fill.1) != 2] <- "."
  tail(sort(table(fill.1)))

  sites = IN.seq$pos
  editseq = fill.1
  sites.temp <- c(0, sites, max(sites)+1)

  for (j in 1:(length(sites.temp)-1)) {
    tofill <- seq(sites.temp[j]+1,(sites.temp[j+1]-1))
    s1 <- editseq[pmax(1, sites.temp[j])]
    s2 <- editseq[pmin(length(editseq), sites.temp[j+1])]

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