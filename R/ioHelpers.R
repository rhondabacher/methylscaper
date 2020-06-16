#' Reading methylation data files
#'
#' @param filepath A string of the file path.
#' @param ... Additional parameters passed to \code{read.table}.
#' @importFrom utils read.table
#' @export
readMethylationData <- function(filepath, ...)
{
  if (grepl(pattern = ".tsv", x = filepath, fixed = TRUE) | 
      grepl(pattern = ".txt", x = filepath, fixed=TRUE))
    read.table(filepath, header=T, row.names = 1,
              stringsAsFactors = F, quote = "", sep = "\t", 
              comment.char = "", ...)
  else
    read.table(filepath, header=T, row.names = 1,
               stringsAsFactors = F, quote = "", sep = ",", 
               comment.char = "", ...)
}

#' Writing methylation data matrices
#' 
#' @param dat A methylation data matrix, output from the \code{runAlign} function.
#' @param filepath A string of the file path.
#' @param ... Additional parameters passed to \code{write.table}.
#' @importFrom utils write.table
#' @export
writeMethylationData <- function(dat, filepath, ...)
{
  if (grepl(pattern = ".csv", x = filepath, fixed=TRUE))
    write.table(dat, file=filepath, quote=F, row.names = F, sep=",", ...)
  else if (grepl(pattern = ".tsv", x = filepath, fixed=TRUE) | 
           grepl(pattern = ".txt", x = filepath, fixed=TRUE))
    write.table(dat, file=filepath, quote=F, row.names = F, sep="\t", ...)
  else # if the file isn't csv, tsv, or txt, we force it to be csv
  {
    filepath <- paste0(filepath, ".csv")
    write.table(align.out$hcg, file=filepath, quote=F, row.names = F, sep=",", ...)
  }
}