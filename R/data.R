#' Example single-molecule GCH and HCG matrices
#' 
#' This was made after running methylscaper::runAlign(fasta=reads_day7, 
#' ref=ref_seq)
#' @docType data
#' @usage data(day7)
"day7"

#' Example reads from single-molecule experiment
#' 
#' This dataset was loaded into R using seqinr::read.fasta
#' @docType data
#' @usage data(reads_day7)
"reads_day7"

#' Example reference sequence to align reads to from a 
#' single-molecule experiment
#' 
#' This dataset was loaded into R using seqinr::read.fasta
#' @docType data
#' @usage data(ref_seq)
"ref_seq"

#' Example preprocessed single-cell experiment subset
#' 
#' This data is from GSE109262, and has been pre-processed by methylscaper
#' It contains a small subaet of chromosome 19 region from 8956841bp - 8976841bp.
#' @docType data
#' @usage data(singlecell_subset)
"singlecell_subset"


#' Human gene symbols and positions
#'
#' library(biomaRt) #v2.44.4
#' ensembl <- useMart("ensembl")
#' ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#' my_chr <- c(1:22, 'M', 'X', 'Y')
#' hum_bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
#'             filters = 'chromosome_name',
#'             values = my_chr,
#'             mart=ensembl)
#' 
#' @docType data
#' @usage data(hum_bm)
"hum_bm"


#' Mouse gene symbols and positions
#'
#' library(biomaRt) #v2.44.4
#' ensembl <- useMart("ensembl")
#' ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
#' my_chr <- c(1:19, 'M', 'X', 'Y')
#' mouse_bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'mgi_symbol'),
#'             filters = 'chromosome_name',
#'             values = my_chr,
#'             mart=ensembl)
#' 
#' @docType data
#' @usage data(mouse_bm)
"mouse_bm"

