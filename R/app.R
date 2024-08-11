# Run the application
#' methylscaper
#'
#' Runs the methylscaper Shiny app.
#'
#' @return This starts up the shiny app interface for methylscaper.
#' @import shiny
#' @import seriation
#' @importFrom seqinr read.fasta
#' @importFrom shinyjs toggleElement showElement
#' @importFrom shinyFiles shinyDirChoose shinyDirButton getVolumes
#' @importFrom grDevices dev.off pdf png cairo_pdf
#' @importFrom utils write.csv capture.output data
#' @importFrom data.table fread
#' @import seriation
#' @export
#' @examples
#'
#' # methylscaper()
methylscaper <- function(port = 19903, is_local_development = TRUE) {
    # default port is set to 19903
    # if is_local_development is set to false, then a browser window will
    # not open up. We should use FALSE in production server
    options_args <- list(
        shiny.maxRequestSize = 10000 * 1024^5,
        shiny.launch.browser = is_local_development,
        shiny.port = port,
        test.mode = getOption("shiny.testmode", FALSE)
    )
    if (is_local_development == FALSE) {
        options_args$shiny.host <- "0.0.0.0"
    }
    do.call(options, options_args)

    shinyApp(ui = ui, server = server)
}
