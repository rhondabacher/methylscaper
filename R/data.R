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

#' Human gene symbols and positions
#'
#' library(biomaRt) #v2.44.4
#' ensembl <- useMart("ensembl")
#' ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#' my_chr <- c(1:22, 'M', 'X', 'Y')
#' hum.bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
#'             filters = 'chromosome_name',
#'             values = my_chr,
#'             mart=ensembl)
#' 
#' @docType data
#' @usage data(hum.bm)
"hum.bm"


#' Mouse gene symbols and positions
#'
#' library(biomaRt) #v2.44.4
#' ensembl <- useMart("ensembl")
#' ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
#' my_chr <- c(1:19, 'M', 'X', 'Y')
#' mouse.bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'mgi_symbol'),
#'             filters = 'chromosome_name',
#'             values = my_chr,
#'             mart=ensembl)
#' 
#' @docType data
#' @usage data(mouse.bm)
"mouse.bm"

