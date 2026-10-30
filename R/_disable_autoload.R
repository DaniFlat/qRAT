# This file tells shiny::runApp() not to source the files in R/ automatically.
# app.R loads the package with pkgload::load_all(); without this file every R/ file was read twice
# at startup (about 1 s extra and the warning "this directory appears to contain an R package").
