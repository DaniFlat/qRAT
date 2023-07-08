#' @name qRAT - qPCR Relative Expression Analysis Tool
#' @author Daniel Flatschacher
#' @title Run the shiny application#
#' @export qRAT

qRAT <- function(...) {

  shiny::addResourcePath('www',
                         system.file('www',
                                     package = 'qRAT'))
shinyApp(ui = ui, server = server, ...)

}
