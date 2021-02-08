# Run the application
#' methylscaper
#'
#' Runs the methylscaper Shiny app.
#'
#' @return This starts up the shiny app interface for methylscaper.
#' @import shiny
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom utils write.csv
#' @importFrom svglite svglite
#' @export
#' @examples 
#' 
#' # methylscaper()
methylscaper <- function() {
    options(shiny.maxRequestSize = 10000*1024^5)

    shinyApp(ui = ui, server = server)}
