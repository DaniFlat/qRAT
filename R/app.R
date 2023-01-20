#' @name qRAT - qPCR Relative Expression Analysis Tool
#' @author Daniel Flatschacher
#' @title Run the shiny application

qRAT <- function(...) {

addResourcePath(prefix = "www", directoryPath = "./www")
shinyApp(ui = ui, server = server, ...)

}
