#' Example single-molecule GCH and HCG matrices
#' 
#' This was made after running methylscaper::runAlign(fasta=reads.day7, 
#' ref=ref.seq)
#' @docType data
#' @usage data(day7)
"day7"

#' Example reads from single-molecule experiment
#' 
#' This dataset was loaded into R using seqinr::read.fasta
#' @docType data
#' @usage data(reads.day7)
"reads.day7"

#' Example reference sequence to align reads to from a 
#' single-molecule experiment
#' 
#' This dataset was loaded into R using seqinr::read.fasta
#' @docType data
#' @usage data(ref.seq)
"ref.seq"

#' Example GCH file for single-cell experiment
#' 
#' This data is from GSE109262, and has been pre-processed by methylscaper
#' @docType data
#' @usage data(chr19_example_GCH)
"chr19_example_GCH"

#' Example HCG file for single-cell experiment
#' 
#' This data GSE109262, and has been pre-processed by methylscaper
#' @docType data
#' @usage data(chr19_example_HCG)
"chr19_example_HCG"