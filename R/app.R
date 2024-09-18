# Run the application
#' methylscaper
#'
#' Runs the methylscaper Shiny app.
#'
#' @return This starts up the shiny app interface for methylscaper.
#' @import shiny
#' @import shinyFiles
#' @import seriation
#' @importFrom shinyjs disabled useShinyjs enable showElement toggleElement
#' @importFrom seqinr read.fasta
#' @importFrom grDevices dev.off pdf png cairo_pdf
#' @importFrom utils write.csv capture.output data
#' @importFrom data.table fread
#' @export
#' @examples
#'
#' # methylscaper()
methylscaper <- function() {
    options(shiny.maxRequestSize = 10000 * 1024^5)

    shinyApp(ui = ui, server = server)
}
