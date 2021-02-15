#' Reading methylation data files
#'
#' @param filepath A string of the file path.
#' @param ... Additional parameters passed to \code{read.table}.
#' @return This function reads in the GCH and HCG matrices if
#'  they have been previously saved within the Shiny app, or
#'  via the writeMethylationData function. 
#' @importFrom utils read.table
#' @export
#' @examples 
#' 
#' # my.gchdata <- readMethylationData("GCH.csv")
#' # my.hcgdata <- readMethylationData("HCG.csv")
readMethylationData <- function(filepath, ...)
{
    if (grepl(pattern = ".tsv", x = filepath, fixed = TRUE) |
        grepl(pattern = ".txt", x = filepath, fixed=TRUE))
    read.table(filepath, header=TRUE, row.names = 1,
                stringsAsFactors = FALSE, quote = "", sep = "\t",
                comment.char = "", ...)
    else
        read.table(filepath, header=TRUE, row.names = 1,
                    stringsAsFactors = FALSE, quote = "", sep = ",",
                    comment.char = "", ...)
}

#' Writing methylation data matrices
#'
#' @param dat A methylation data matrix, output from the
#'  \code{runAlign} function.
#' @param filepath A string of the file path.
#' @param ... Additional parameters passed to \code{write.table}.
#' @return This function writes the GCH and HCG matrices as a csv file
#'  by default, although a tab separated file is also valid and can be 
#'  specified as .txt or .tsv.
#' @importFrom utils write.table
#' @examples 
#' 
#' data(day7)
#' # writeMethylationData(day7$gch, "GCH.csv")
#' 
#' @export
writeMethylationData <- function(dat, filepath, ...)
{
    if (grepl(pattern = ".csv", x = filepath, fixed=TRUE))
        write.table(dat, file=filepath, quote=FALSE, 
                row.names = FALSE, sep=",", ...)
    else if (grepl(pattern = ".tsv", x = filepath, fixed=TRUE) |
            grepl(pattern = ".txt", x = filepath, fixed=TRUE))
    write.table(dat, file=filepath, quote=FALSE, 
            row.names = FALSE, sep="\t", ...)
    else # if the file isn't csv, tsv, or txt, we force it to be csv
    {
    filepath <- paste0(filepath, ".csv")
    write.table(dat, file=filepath, quote=FALSE, 
            row.names = FALSE, sep=",", ...)
    }
}
