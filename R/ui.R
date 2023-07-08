#' @name qRAT - qPCR Relative Expression Analysis Tool
#' @title UI for the qRAT Shiny app
#' @author Daniel Flatschacher
####

####
# Packages
####

library("shiny")
library("bslib")
library("shinyWidgets")
library("HTqPCR")
library("ddCt")
library("plotly")
library("shinyjs")
library("ggplot2")
library("scales")
library("data.table")
library("DT")
library("dplyr")
library("waiter")
library("tidyr")
library("stringr")
library("magrittr")
library("shinycssloaders")
library("viridisLite")

####
# Loading Screen
####

options(spinner.color = "#29abe0", spinner.type = "8")

# color list for plots
colorlist <- c("Paired", "Spectral", "Blues", "Greens", "Greys", "Oranges", "Purples", "Reds", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "viridis")

#shinyLink functions
shinyLink <- function(to, label) {
  tags$a(
    class = "shiny__link",
    href = to,
    label,
    class = "btn btn-success",
    role = "button"
  )
}

shinyLinkInText <- function(to, label) {
  tags$a(
    class = "shiny__link",
    href = to,
    label
  )
}



####
# User Interface
####
# Define UI
ui <- page_navbar(id = "navbar",

  # setup of theme
  theme = bs_theme(
    bootswatch = "flatly",
    version = 5,
    font_scale = 1,
    base_font = "Helvetica",
    success = "#93c54b",
    warning = "#f47c3c",
    info = "#29abe0",
    primary = "#325d88",
  ),
  title = "qRAT",
  collapsible = TRUE,
  inverse = FALSE,
  fillable = FALSE,
  tabPanel("Start",
    icon = icon("house"), class = "active",

    # loading screen
    waiter::use_waiter(),
    waiter::use_hostess(),
    waiter::use_steward(colors = c("#325d88", "#2c3e50", "#325d88", "#6610f2")),
    waiter::waiter_preloader(html = tagList(waiter::spin_flower(), HTML(paste(tags$span(style = "font-size:17px;letter-spacing: 1.5px", "Loading"), tags$span(style = "font-size:17px;font-weight:600", "qRAT..."), sep = " ")))),
    tags$script(src = "www/shinyLink.js"), # shiny Link for internal linking (to tabsets)
    tags$style(HTML(".help-block {color: #325d88 !important;}")), # color of helpText

    # Start Page HTML
    div(
      class = "jumbotron",
      div(
        class = "container",
        div(
          class = "row",
          div(
            class = "col-md-8",
            h1(class = "display-3", "qRAT"),
            p(class = "lead", "qPCR - Relative Expression Analysis Tool"),
            hr(class = "my-4", "qRAT is an R based Electron app designed to provide a comprehensive, simple and easy to use tool for real-time quantitative PCR data analysis.
     Functions are provided for users who have little statistical and no R programming background."),
            br(),
            br(),
            p(
              class = "lead",
              a(class = "btn btn-primary btn-lg", href = "https://uibk.ac.at/microbiology/services/qrat/", target = "_blank", role = "button", "Learn more")
            )
          ),
          div(
            class = "col-md-4",
            img(src = "www/logo.svg", class = "img-fluid", width = "80%", height = "80%")
          )
        )
      )
    ),
    p(),
    fluidRow(
      column(
        4,
        div(
          class = "card text-white bg-success mb-3",
          div(class = "card-header", "Directions"),
          div(
            class = "card-body",
            h4(class = "card-title", "How To?"),
            p(class = "card-text", "- Extensive documentation on data input, manipulation and analysis"),
            p(class = "card-text", "- Checklist and most important points")
          )
        ),
        br(),
        p(
          class = "lead",
          shinyLink(to = "help", label = "Learn more"),
        )
      ),
      column(
        4,
        div(
          class = "card text-white bg-warning mb-3",
          div(class = "card-header", "App Status"),
          div(
            class = "card-body",
            h4(class = "card-title", "qRAT version"),
            p(class = "card-text", "You're running: Version 0.2.0"),
            p(class = "card-text", "It is recommended to check for Updates before using the application")
          )
        ),
        br(),
        p(
          class = "lead",
          actionButton(class = "btn btn-warning", label = "Check for Update!", inputId = "updateCheck")
        )
      ),
      column(
        4,
        div(
          class = "card text-white bg-info mb-3",
          div(class = "card-header", "Citation"),
          div(
            class = "card-body",
            h4(class = "card-title", "Note"),
            p(class = "card-text", "If you use qRAT in published research, please include the correct citations."),
            p(class = "card-text", "Thanks!")
          )
        ),
        br(),
        p(
          class = "lead",
          actionButton(class = "btn btn-info", label = "How to cite?", inputId = "Citation")
        )
      ),
    ),
    p(),
    p(),
    p(class = "text-muted text-center", "Copyright ©2021-2023 Daniel Flatschacher, Department of Microbiology, University of Innsbruck"),
  ),
  tabPanel("Single Plate",
    icon = icon("stop"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        conditionalPanel(
          condition = "input.tabsSingle=='Raw Data'",
          h4("Single Plate Input"),
          helpText("Choose the separator used in the file. If you are not sure, set it to 'auto'.", tags$br(),
                   "See", shinyLinkInText(to = "help", label = "User Guide"), "for more details."),
          radioGroupButtons("Sep",
            label = "Data separator", choices = c("auto", "comma", "tab", "semicolon"),
            checkIcon = list(
              yes = icon("ok",
                lib = "glyphicon"
              ),
              no = icon("remove",
                lib = "glyphicon"
              )
            )
          ),
          helpText("Choose your file", tags$b("(csv or txt)"), tags$br(),
                   "See", shinyLinkInText(to = "help", label = "Help"), "to download example data file."),
          fileInput("dtfile", label = "Upload single plate data file here", accept = c(
            "text/csv", "text/comma-separated-values", "text/tab-separated-values",
            "text/plain", ".csv", ".txt"
          )),
        ),
        conditionalPanel(
          condition = "input.tabsSingle=='Filtering & Quality'",
          h4("Analysis Input"),
          h5("Filtering"),
          h5("Define Outlier (optional)"),
          actionLink(icon = icon("circle-info"), label=NULL, style="color: #325d88", inputId = "outlierSP"),
          helpText("Set a range of acceptable Cq values"),
          sliderInput("maxCt", "Cq Cut-off", min = 0, max = 45, value = c(5, 35), step = 1),
          helpText("Set max deviation between replicates (from mean)"),
          sliderInput("repDev", "Replicate Variability", min = 0, max = 2, value = 0.3, step = 0.01),
          h5("Experimental Controls"),
          helpText("Choose NTC(s) or RT- to be excluded from calculations. If the sample set does not include NTC(s) or RT- leave Input blank."),
          virtualSelectInput(inputId = "NTC_Input",
                      label = "Select NTC and RT- Samples",
                      choices = attr("NTCs", "Labels"),
                      multiple = TRUE),
        ),
        conditionalPanel(
          condition = "input.tabsSingle=='Relative dCq' || input.tabsSingle=='Relative ddCq'",
          h4("Analysis Input"),
          h5("Housekeeping Gene(s)", actionLink(icon = icon("circle-info"), label=NULL, style="color: #325d88", inputId = "housekeepingInfo")),
          helpText("Set reference gene(s)"),
          virtualSelectInput("Refs",
            label = "Reference Genes",
            choices = attr("Refs", "Labels"),
            multiple = TRUE
          )
        ),
        conditionalPanel(
          condition = "input.tabsSingle=='Relative dCq'",
          h5("Plot Settings"),
          helpText("Show Settings"),
          switchInput(inputId = "show_plot_settings_dCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_plot_settings_dCq == true",
            virtualSelectInput(inputId = "PlotData", label = "Plot Data", choices = c("RQ", "dCq", "-dCq")),
            virtualSelectInput(inputId = "SamplePicker", label = "Samples", choices = "", multiple = TRUE),
            virtualSelectInput(inputId = "GenePicker", label = "Genes", choices = "", multiple = TRUE),
            virtualSelectInput(inputId = "PlotType", label = "Plot Type", choices = c("Bar Chart", "Dot Plot")),
            checkboxInput("ShowErrorBar", "Show error bars", FALSE),
            virtualSelectInput(inputId = "scale", label = "Scale", choices = c("normal", "log"))),

          h5("Customize Plot Appearance"),
          helpText("Show Appearance Settings"),
          switchInput(inputId = "show_appearance_settings_dCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_appearance_settings_dCq == true",
            virtualSelectInput(inputId = "colorpicker", label = "Pick Colors", choices = colorlist),
            h6("Define Sample Order"),
            helpText("Drag and Drop Samples"),
            orderInput("x_order", "", items = NULL, width = "100%"),
            h6("X Axis Text Angle"),
            helpText("Rotation angle of x axis labels"),
            sliderTextInput(inputId = "xTextAngle_dCq", label = "", choices = c(0, 30, 45, 60, 90), grid = TRUE),
            textInput(inputId = "plotTitle_dCq", label = "Plot Title", value = "", width = NULL, placeholder = NULL),
            textInput(inputId = "legendTitle_dCq", label = "Legend Title", value = "Genes", width = NULL, placeholder = "Genes"),
            virtualSelectInput(inputId = "legendPosition_dCq", label = "Legend Position", choices = c("right", "none"))),

          h5("Plot Export"),
          helpText("Show Plot Export Settings:"),
          switchInput(inputId = "show_export_settings_dCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_export_settings_dCq == true",
            virtualSelectInput(inputId = "exportFormat", label = "File Format", choices = c("svg", "png", "jpeg", "webp")),
            numericInput("width_dCq", label = "Width", value = 800),
            numericInput("height_dCq", label = "Height", value = 800),
            numericInput("scale_dCq", label = "Scale", value = 3),
            helpText("set scale to 6 for 600 dpi"))
        ),
        conditionalPanel(
          condition = "input.tabsSingle=='Relative ddCq'",
          h5("Calibrator"),
          helpText("Set Calibrator for ddCq"),
          virtualSelectInput("Mock", label = "Calibrator", choices = ""),
          h5("Plot Settings"),
          helpText("Show Settings"),
          switchInput(inputId = "show_plot_settings_ddCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_plot_settings_ddCq == true",
            virtualSelectInput(inputId = "PlotDataDDCt", label = "Plot Data", choices = c("Fold Change", "ddCt")),
            virtualSelectInput("SamplePickerDDCt", label = "Samples", choices = "", multiple = TRUE),
            virtualSelectInput("GenePickerDDCt", label = "Genes", choices = "", multiple = TRUE),
            virtualSelectInput(inputId = "PlotTypeDDCt", label = "Plot Type", choices = c("Bar Chart", "Dot Plot")),
            checkboxInput("ShowErrorBarDDCt", "Show error bars", FALSE),
            virtualSelectInput(inputId = "scaleDDCt", label = "Scale", choices = c("normal", "log"))),

          h5("Customize Plot Appearance"),
          helpText("Show Appearance Settings"),
          switchInput(inputId = "show_appearance_settings_ddCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_appearance_settings_ddCq == true",
            virtualSelectInput(inputId = "colorpickerDDCt", label = "Pick Colors", choices = colorlist),
            h6("Define Sample Order"),
            helpText("Drag and Drop Samples"),
            orderInput("x_orderDDCt", "", items = NULL, width = "100%"),
            h6("X Axis Text Angle"),
            helpText("Rotation angle of x axis labels"),
            sliderTextInput(inputId = "xTextAngle_ddCq", label = "", choices = c(0, 30, 45, 60, 90), grid = TRUE),
            textInput(inputId = "plotTitle_ddCq", label = "Plot Title", value = "", width = NULL, placeholder = NULL),
            textInput(inputId = "legendTitle_ddCq", label = "Legend Title", value = "Genes", width = NULL, placeholder = "Genes"),
            virtualSelectInput(inputId = "legendPosition_ddCq", label = "Legend Position", choices = c("right", "none"))),

          h5("Plot Export"),
          helpText("Show Plot Export Settings:"),
          switchInput(inputId = "show_export_settings_ddCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_export_settings_ddCq == true",
            virtualSelectInput(inputId = "exportFormatDDCt", label = "File Format", choices = c("svg", "png", "jpeg", "webp")),
            numericInput("width_ddCq", label = "Width", value = 800),
            numericInput("height_ddCq", label = "Height", value = 800),
            numericInput("scale_ddCq", label = "Scale", value = 3),
            helpText("set scale to 6 for 600 dpi"))
        ),
        conditionalPanel(
          condition = "input.tabsSingle=='Statistical Analysis'",
          h4("Analysis Input"),
          virtualSelectInput(inputId = "limma_input", label = "Limma Input", choices = c("dCq")),
          h5("Statistical Analysis"),
          helpText("The sample names must be syntactically valid variable names in R and so, for example, must
          begin with a letter rather than a numeral. Special characters (also '.', '_', etc. should be avoided)"),
          h5("Parameters"),
          helpText("Adjust p-values"),
          virtualSelectInput(inputId = "adjustMethod", label = "Adjustment Method", choices = c("Benjamini & Hochberg", "Holm", "Bonferroni")),
          helpText("Choose your type of comparison"),
          radioButtons("compType", label = "Comparison Type", choices = c("Single Comparison", "Multiple paired Comparisons"), inline = TRUE),

          #### Single Comparison Settings Start
          conditionalPanel(
            condition = "input.compType=='Single Comparison'",
            helpText("Choose a single control:"),
            virtualSelectInput("Comps", label = "Single Comparison", choices = "")
          ),
          #### Single Comparison Settings End

          #### Multiple Comparisons Settings Start
          conditionalPanel(
            condition = "input.compType=='Multiple paired Comparisons'",
            h5("Multiple Comparisons"),
            helpText("Add as many paired comparisons as you like. Choose two samples you want to compare:"), br(),
            # helpText("Number of comparisons:"), textOutput("comparisonBoxes_counter"),
            actionButton("add_btn", "Add Comparison", class = "btn-primary"), br(),
            actionButton("rm_btn", "Remove Comparison", class = "btn-danger"), br(),
            uiOutput("comparisonBoxes_ui")
          )
          #### Multiple Comparisons Settings End
        ),
        conditionalPanel(
          condition = "input.tabsSingle=='MIQE Check'",
          h5("MIQE"),
          helpText("What is MIQE?"),
          p(
            "The aim of the", strong("M"), "inimum", strong("I"), "nformation for Publication of", strong("Q"), "uantitative Real-Time PCR", strong("E"), "xperiments guidelines", a(href = "http://rdml.org/miqe.html", "(MIQE & RDML)"),
            "is to provide authors, reviewers and editors with specifications for the minimum information that must be reported for a qPCR experiment."
          ),
          p("MIQE ensures its relevance, accuracy, correct interpretation and repeatability.")
        ),
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs", id = "tabsSingle",
          tabPanel(
            "Raw Data", h4("Raw Data", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedSingle",
              fluidRow(
                column(5, h5("Data Table"), DTOutput("dataSinglePlate")),
                column(6, h5(""), plotlyOutput("plotCtCard"), plotlyOutput("plotCtDistrib"))
              )
            )
          ),
          tabPanel(
            "Filtering & Quality", h4("Filtering & Quality", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedSingle",
              fluidRow(
                column(12,
                       fluidRow(
                        column(12, h5("Wells to be excluded from analysis"), DTOutput("dataSinglePlateBadRep")), br(), br()
                        )),
              ),
              fluidRow(column(6, h5(""), plotlyOutput("SinglePlateFilterPlot")),
                       column(6, h5(""), uiOutput("SinglePlateMIQEcheck"))),
              fluidRow(column(12, h5(""), plotlyOutput("SinglePlateBoxplot")))
            )
          ),
          tabPanel(
            "Relative dCq", h4("dCq Analysis", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedSingle",
              fluidRow(h5("dCq Graph Output"), plotlyOutput("ddctAbsGraph", width = "100%", height = "50%")), br()
            ),
            h5("dCq Table Output"), fluidRow(DTOutput("ddctAbsolute")),
          ),
          tabPanel(
            "Relative ddCq", h4("ddCq Analysis", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedSingle",
              fluidRow(h5("ddCq Graph Output"), plotlyOutput("ddctRelGraph", width = "100%", height = "50%")), br()
            ),
            h5("ddCq Table Output"), fluidRow(DTOutput("ddctRelative")),
          ),
          #    tabPanel(
          #      "MIQE Check", h4("MIQE Guidline Check", align = "center"), br(), br(),
          #      conditionalPanel(
          #      condition = "output.fileUploadedSingle",
          #      fluidRow(
          #        column(6, uiOutput("SinglePlateMIQEcheck")),
          #        column(
          #          6, p(""),
          #          strong("Appropiate number of technical replicates:"), "n >= 2"
          #        )
          #      ))
          #    ),
          tabPanel(
            "Statistical Analysis", h4("Statistical Analysis", align = "center"), h5("Function for detecting differentially expressed genes"),
            br(),
            DTOutput("resultLimma")
          )
        )
      )
    )
  ),
  tabPanel("Multiple Plates",
    icon = icon("table-cells-large"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        conditionalPanel(
          condition = "input.tabsMulti=='Raw Data'",
          h4("Multiple Plates Input"),
          helpText("Choose the separator used in the file. If you are not sure, set it to 'auto'.", tags$br(),
                   "See", shinyLinkInText(to = "help", label = "User Guide"), "for more details."),
          radioGroupButtons("SepM",
            label = "Data separator", choices = c("auto", "comma", "tab", "semicolon"),
            checkIcon = list(
              yes = icon("ok",
                lib = "glyphicon"
              ),
              no = icon("remove",
                lib = "glyphicon"
              )
            )
          ),
          helpText("Choose your files", tags$b("(csv or txt)"), tags$br(),
                   "Multiple files can be selected and loaded while holding the [Ctrl] key.", tags$br(),
                   "See", shinyLinkInText(to = "help", label = "User Guide"), "for more details."),
          fileInput("plates", label = "Upload multiple plates data files here", accept = c(
            "text/csv", "text/comma-separated-values", "text/tab-separated-values",
            "text/plain", ".csv", ".txt"
          ), multiple = TRUE)
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='Filtering & Quality'",
          h4("Analysis Input"),
          h5("Filtering"),
          helpText("Set a range of acceptable Cq values"),
          sliderInput("maxCtMulti", "Cq Cut-off", min = 0, max = 45, value = c(5, 35), step = 1),
          helpText("Set max deviation between replicates (from mean)"),
          sliderInput("repDevMulti", "Replicate Variability", min = 0, max = 2, value = 0.3, step = 0.01),
          virtualSelectInput(inputId = "NTC_Input_MP",
                      label = "Select NTC and RT- Samples",
                      choices = attr("MP_NTCs", "Labels"),
                      multiple = TRUE),
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='Inter-Plate Calibration'",
          h4("Analysis Input"),
          h5("Inter-Plate Calibrator"),
          helpText("Choose if you want to use Inter-Plate Calibration function"),
          prettySwitch(
            inputId = "Id027",
            label = "Use Inter-Plate Calibration",
            status = "success",
            fill = TRUE
          ),
          helpText("Set Sample as Inter-Plate calibrator"),
          virtualSelectInput("IPC", label = "IPC", choices = ""),
          helpText("Choose Samples for Comparison-Plot"),
          virtualSelectInput("SamplePickerIPCcomparison",
                      label = "Samples to Plot",
                      choices = "",
                      multiple = TRUE
          ),
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='Relative dCq' || input.tabsMulti=='Relative ddCq'",
          h4("Analysis Input"),
          h5("Housekeeping Gene(s)", actionLink(icon = icon("circle-info"), label=NULL, style="color: #325d88", inputId = "housekeepingInfoMP")),
          helpText("Set reference gene(s)"),
          virtualSelectInput("RefsM",
            label = "Reference Genes",
            choices = attr("Refs", "Labels"),
            multiple = TRUE
          ),
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='Relative dCq'",
          h5("Plot Settings"),
          helpText("Show Settings"),
          switchInput(inputId = "show_plot_settings_MultidCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_plot_settings_MultidCq == true",
            virtualSelectInput(inputId = "PlotDataMulti", label = "Plot Data", choices = c("RQ", "dCq", "-dCq")),
            virtualSelectInput("SamplePickerMulti", label = "Samples", choices = "", multiple = TRUE),
            virtualSelectInput("GenePickerMulti", label = "Genes", choices = "", multiple = TRUE),
            virtualSelectInput(inputId = "PlotTypeMulti", label = "Plot Type", choices = c("Bar Chart", "Dot Plot")),
            checkboxInput("ShowErrorBarMulti", "Show error bars", FALSE),
            virtualSelectInput(inputId = "scaleMulti", label = "Scale", choices = c("normal", "log"))),

          h5("Customize Plot Appearance"),
          helpText("Show Appearance Settings"),
          switchInput(inputId = "show_appearance_settings_MultidCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_appearance_settings_MultidCq == true",
            virtualSelectInput(inputId = "colorpickerMulti", label = "Pick Colors", choices = colorlist),
            h6("Define Sample Order"),
            helpText("Drag and Drop Samples"),
            orderInput("x_order_MultidCq", "", items = NULL, width = "100%"),
            h6("X Axis Text Angle"),
            helpText("Rotation angle of x axis labels"),
            sliderTextInput(inputId = "xTextAngle_MultidCq", label = "", choices = c(0, 30, 45, 60, 90), grid = TRUE),
            textInput(inputId = "plotTitle_MultidCq", label = "Plot Title", value = "", width = NULL, placeholder = NULL),
            textInput(inputId = "legendTitle_MultidCq", label = "Legend Title", value = "Genes", width = NULL, placeholder = "Genes"),
            virtualSelectInput(inputId = "legendPosition_MultidCq", label = "Legend Position", choices = c("right", "none"))),

          h5("Plot Export"),
          helpText("Show Plot Export Settings:"),
          switchInput(inputId = "show_export_settings_MultidCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_export_settings_MultidCq == true",
            virtualSelectInput(inputId = "exportFormatMulti", label = "File Format", choices = c("svg", "png", "jpeg", "webp")),
            numericInput("width_dCqMulti", label = "Width", value = 800),
            numericInput("height_dCqMulti", label = "Height", value = 800),
            numericInput("scale_dCqMulti", label = "Scale", value = 3),
            helpText("set scale to 6 for 600 dpi"))
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='Relative ddCq'",
          h5("Calibrator"),
          helpText("Set calibrator for ddCq"),
          virtualSelectInput("MockM", label = "Calibrator", choices = ""),
          h5("Plot Settings"),
          helpText("Show Settings"),
          switchInput(inputId = "show_plot_settings_MultiddCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_plot_settings_MultiddCq == true",
            virtualSelectInput(inputId = "PlotDataDDCtMulti", label = "Plot Data", choices = c("Fold Change", "ddCt")),
            virtualSelectInput("SamplePickerDDCtMulti", label = "Samples", choices = "", multiple = TRUE),
            virtualSelectInput("GenePickerDDCtMulti", label = "Genes", choices = "", multiple = TRUE),
            virtualSelectInput(inputId = "PlotTypeDDCtMulti", label = "Plot Type", choices = c("Bar Chart", "Dot Plot")),
            checkboxInput("ShowErrorBarDDCtMulti", "Show error bars", FALSE),
            virtualSelectInput(inputId = "scaleDDCtMulti", label = "Scale", choices = c("normal", "log"))),

          h5("Customize Plot Appearance"),
          helpText("Show Appearance Settings"),
          switchInput(inputId = "show_appearance_settings_MultiddCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_appearance_settings_MultiddCq == true",
            virtualSelectInput(inputId = "colorpickerDDCtMulti", label = "Pick Colors", choices = colorlist),
            h6("Define Sample Order"),
            helpText("Drag and Drop Samples"),
            orderInput("x_order_MultiddCq", "", items = NULL, width = "100%"),
            h6("X Axis Text Angle"),
            helpText("Rotation angle of x axis labels"),
            sliderTextInput(inputId = "xTextAngle_MultiddCq", label = "", choices = c(0, 30, 45, 60, 90), grid = TRUE),
            textInput(inputId = "plotTitle_MultiddCq", label = "Plot Title", value = "", width = NULL, placeholder = NULL),
            textInput(inputId = "legendTitle_MultiddCq", label = "Legend Title", value = "Genes", width = NULL, placeholder = "Genes"),
            virtualSelectInput(inputId = "legendPosition_MultiddCq", label = "Legend Position", choices = c("right", "none"))),

          h5("Plot Export"),
          helpText("Show Plot Export Settings:"),
          switchInput(inputId = "show_export_settings_MultiddCq",onLabel = "Yes", offLabel = "No"),
          conditionalPanel(
            condition = "input.show_export_settings_MultiddCq == true",
            virtualSelectInput(inputId = "exportFormatDDCtMulti", label = "File Format", choices = c("svg", "png", "jpeg", "webp")),
            numericInput("width_ddCqMulti", label = "Width", value = 800),
            numericInput("height_ddCqMulti", label = "Height", value = 800),
            numericInput("scale_ddCqMulti", label = "Scale", value = 3),
            helpText("set scale to 6 for 600 dpi"))
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='Statistical Analysis'",
          h4("Analysis Input"),
          virtualSelectInput(inputId = "limma_inputMulti", label = "Limma Input", choices = c("dCq")),
          h5("Statistical Analysis"),
          helpText("The sample names must be syntactically valid variable names in R and so, for example, must
          begin with a letter rather than a numeral. Special characters (also '.', '_', etc. should be avoided)"),
          h5("Parameters"),
          helpText("Adjust p-values"),
          virtualSelectInput(inputId = "adjustMethodMulti", label = "Adjustment Method", choices = c("Benjamini & Hochberg", "Holm", "Bonferroni")),
          helpText("Choose your type of comparison"),
          radioButtons("compTypeM", label = "Comparison Type", choices = c("Single Comparison", "Multiple paired Comparisons"), inline = TRUE),

          #### Single Comparison Settings Start
          conditionalPanel(
            condition = "input.compTypeM=='Single Comparison'",
            helpText("Choose a single control:"),
            virtualSelectInput("CompsM", label = "Single Comparison", choices = "")
          ),
          #### Single Comparison Settings End

          #### Multiple Comparisons Settings Start
          conditionalPanel(
            condition = "input.compTypeM=='Multiple paired Comparisons'",
            h5("Multiple Comparisons"),
            helpText("Add as many paired comparisons as you like. Choose two samples you want to compare:"), br(),
            actionButton("add_btnM", "Add Comparison", class = "btn-primary"), br(),
            actionButton("rm_btnM", "Remove Comparison", class = "btn-danger"), br(),
            uiOutput("comparisonBoxesM_ui")
          )
          #### Multiple Comparisons Settings End
        ),
        conditionalPanel(
          condition = "input.tabsMulti=='MIQE Check'",
          h5("MIQE"),
          helpText("What is MIQE?"),
          p(
            "The aim of the", strong("M"), "inimum", strong("I"), "nformation for Publication of", strong("Q"), "uantitative Real-Time PCR", strong("E"), "xperiments guidelines", a(href = "http://rdml.org/miqe.html", "(MIQE & RDML)"),
            "is to provide authors, reviewers and editors with specifications for the minimum information that must be reported for a qPCR experiment."
          ),
          p("MIQE ensures its relevance, accuracy, correct interpretation and repeatability.")
        ),
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs", id = "tabsMulti",
          tabPanel(
            "Raw Data", h4("Raw Data", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedMulti",
              fluidRow(
                column(6, h5("Data Table"), DTOutput("multiplePlatesData") %>% withSpinner()),
                column(6, h5(""), align="center", plotlyOutput("MultiplotCtCard"), virtualSelectInput("PlateSelect", label = "Choose Plate", choices = "", width = "20%"), plotlyOutput("MultiplotCtDistrib"))
              )
            )
          ),
          tabPanel(
            "Filtering & Quality", h4("Filtering & Quality", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedMulti",
              fluidRow(
                column(12, h5("Wells to be excluded from analysis"), DTOutput("dataMultiplePlatesBadRep"))), br(), br(),
              fluidRow(
                column(6, h5(""), plotlyOutput("MultiplePlatesFilterPlot")),
                column(6, h5(""), uiOutput("MultiplePlateMIQEcheck"))), br(),
              fluidRow(column(12, h5(""), plotlyOutput("MultiplePlatesBoxplot")))
            )
          ),
          tabPanel(
            "Inter-Plate Calibration", h4("Inter-Plate Calibration", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedMulti",
              fluidRow(column(12, h5("Calibrated Plates"), DTOutput("interPlateCalibration")),
              fluidRow(column(8, h5("Inter-Plate Calibrators"), DTOutput("extractedIPCsTable")),
                       column(4, h5("Calibration Factors"), DTOutput("tableCalibrationFactors")),
              )
              ), br(), br(),
              fluidRow(column(12, plotlyOutput("comparisonCalibrated")))
            )
          ),
          tabPanel(
            "Relative dCq", h4("Relative dCq", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedMulti",
              fluidRow(h5("dCq Graph Output"), plotlyOutput("ddctMultiAbsGraph", width = "100%", height = "50%")), br(),
              h5("dCq Table Output"), fluidRow(DTOutput("ddctAbsoluteMulti"))
            )
          ),
          tabPanel(
            "Relative ddCq", h4("Relative ddCq", align = "center"),
            conditionalPanel(
              condition = "output.fileUploadedMulti",
              fluidRow(h5("ddCq Graph Output"), plotlyOutput("ddctRelGraphMulti", width = "100%", height = "50%")), br(),
              h5("ddCq Table Output"), DTOutput("ddctRelativeMulti")
            )
          ),
          #     tabPanel(
          #       "MIQE Check", h4("MIQE Guidline Check", align = "center"), br(), br(),
          #       conditionalPanel(
          #          condition = "output.fileUploadedMulti",
          #         fluidRow(
          #            column(6,),
          #           column(
          #             6, p(""),
          #              strong("Appropiate number of technical replicates:"), "n >= 2"
          #           )
          #          ))
          #     ),
          tabPanel(
            "Statistical Analysis", h4("Statistical Analysis", align = "center"), h5("Function for detecting differentially expressed genes"),
            br(),
            DTOutput("resultLimmaMulti")
          )
        )
      )
    )
  ),
  tabPanel("Help & About",
    icon = icon("circle-question"),
    value = "help",
    includeHTML("help.html")
  )
)
