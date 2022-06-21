# Run the application 
qRAT <- function(...) {

addResourcePath(prefix = "www", directoryPath = "./www")
shinyApp(ui = ui, server = server, ...)

}