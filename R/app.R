# Run the application
#' methylScaper
#'
#' Runs the methylScaper Shiny app.
#'
#' @import shiny
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom utils write.csv
#' @importFrom svglite svglite
#' @export
methylScaper <- function() {

	options(shiny.maxRequestSize = 10000*1024^5) 

	shinyApp(ui = ui, server = server)}
