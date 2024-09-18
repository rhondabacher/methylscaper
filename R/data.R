#' Example preprocessed single-molecule experiment
#'
#' The RDS in ext data was made specifically with the command:
#' singlemolecule_example <- methylscaper::runAlign(fasta=reads_sm,
#' ref=ref_seq)
#' saveRDS(singlemolecule_example, file="methylscaper/inst/ext/singlemolecule_example.rds", compress = 'xz')
#' A version is also saved as RData used running examples in the man pages.
#' save(singlemolecule_example, file="methylscaper/data/singlemolecule_example.RData", compress = 'xz')
#' @docType data
#' @usage data(singlemolecule_example)
"singlemolecule_example"

#' Example reads from single-molecule experiment
#'
#' This dataset was loaded into R using seqinr::read.fasta
#' @docType data
#' @usage data(reads_sm)
"reads_sm"

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
#' It contains a small subset of chromosome 19 region from 8947041bp - 8987041bp.
#' The RDS in ext data was made specifically with the two commands:
#' singlecell_subset <- subsetSC("~/Downloads/GSE109262_RAW/", chromosome="19",
#'    startPos = 8967041-20000, endPos = 8967041+20000, updateProgress = NULL)
#' saveRDS(singlecell_subset, file="methylscaper/inst/ext/singlecell_subset.rds", compress = 'xz')
#' A version is also saved as RData used running examples in the man pages.
#' save(singlecell_subset, file="methylscaper/data/singlecell_subset.RData", compress = 'xz')
#' @docType data
#' @usage data(singlecell_subset)
"singlecell_subset"

#' Human gene symbols and positions
#'
#' library(biomaRt) #v2.44.4
#' ensembl <- useMart("ensembl") # GRCh38
#' ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#' my_chr <- c(1:22, 'M', 'X', 'Y')
#' human_bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
#'             filters = 'chromosome_name',
#'             values = my_chr,
#'             mart=ensembl)
#'
#' @docType data
#' @usage data(human_bm)
"human_bm"

#' Mouse gene symbols and positions
#'
#' library(biomaRt) #v2.44.4
#' ensembl <- useMart("ensembl") # GRCm39
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
