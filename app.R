source("R/ui.R")
source("R/server.R")
source("R/shinyLink.R")
source("R/qPCR.R")
source("R/PlotSettings.R")

shinyApp(ui = ui, server = server)
