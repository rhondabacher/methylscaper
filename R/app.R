# Run the application
#' methylscaper
#'
#' Runs the methylscaper Shiny app.
#'
#' @import shiny
#' @importFrom shinyFiles shinyDirChoose shinyDirButton getVolumes
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom utils write.csv
#' @importFrom svglite svglite
#' @importFrom data.table fread
#' @importFrom Rcpp sourceCpp
#' @useDynLib methylscaper, .registration = TRUE
#' @export
methylscaper <- function() {

	options(shiny.maxRequestSize = 10000*1024^5)

	shinyApp(ui = ui, server = server)}
