#' @name qRAT - qPCR Relative Expression Analysis Tool
#' @title UI for the qRAT Shiny app
#' @author Daniel Flatschacher

library("shiny")
library("bslib")
library("shinyWidgets")
library("ggiraph")
library("shinyjs")
library("scales")
library("DT")
library("waiter")
library("shinycssloaders")
library("viridisLite")

####
# Configuration & Design Tokens
####
accent_blue <- "#325d88"
primary_col <- "#5a7d9a" 
success_col <- "#4e9a5a"
warning_col <- "#f47c3c"
info_col    <- "#29abe0"
page_bg_col  <- "#ffffff" 
card_header_bg <- "#ecf0f1" 
card_bg_col  <- "#ffffff" 

options(spinner.color = primary_col, spinner.type = "8")
colorlist <- c("Paired", "Spectral", "Blues", "Greens", "Greys", "Oranges", "Purples", "Reds", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdYlBu", "RdYlGn", "viridis")


####
# Bar below the Plate View: wells clicked in the heatmap are applied as manual exclusions.
####
plate_selection_bar <- function(button_id, note = NULL) {
  div(
    class = "d-flex flex-wrap align-items-center gap-2 mt-2",
    actionButton(button_id, "Exclude selected wells", icon = icon("ban"), class = "btn-outline-danger btn-sm"),
    span(class = "text-muted small",
         "Click wells to select them (outlined). Applying switches to Manual Filtering; an empty selection removes the exclusions.",
         note)
  )
}


####
# Sidebar settings of the four plot tabs (dCq/ddCq, Single/Multi Plate) as an accordion.
# Replaces the former "Appearance"/"Export Opt." switches; all input IDs are unchanged.
#   sfx:            suffix of the plot inputs ("", "DDCt", "Multi", "DDCtMulti")
#   exp:            suffix of the export inputs ("dCq", "ddCq", "Multi", "DDCtMulti")
#   download_id:    ID of the download button
#   data_choices:   choices for the data type (RQ/dCq or Fold Change/ddCt)
#   format_choices: export formats (first entry is the default)
#   calibrator:     optional calibrator input (ddCq tabs only)
####
plot_settings_sidebar <- function(sfx, exp, download_id, data_choices, format_choices, calibrator = NULL) {
  accordion(
    multiple = TRUE,
    open = c("Data Selection", "Data Visualization"),
    class = "mb-3",
    accordion_panel(
      "Data Selection", icon = bsicons::bs_icon("funnel"),
      calibrator,
      virtualSelectInput(paste0("SamplePicker", sfx), "Filter Samples", choices = NULL, multiple = TRUE, search = TRUE),
      virtualSelectInput(paste0("GenePicker", sfx), "Filter Genes", choices = NULL, multiple = TRUE, search = TRUE)
    ),
    accordion_panel(
      "Data Visualization", icon = bsicons::bs_icon("bar-chart"),
      layout_column_wrap(
        width = 1/2,
        virtualSelectInput(paste0("PlotData", sfx), "Data Type", choices = data_choices),
        virtualSelectInput(paste0("PlotType", sfx), "Plot Type", choices = c("Bar Chart", "Dot Plot"))
      ),
      layout_column_wrap(
        width = 1/2,
        virtualSelectInput(paste0("scale", sfx), "Y-Axis", choices = c("normal", "log10", "log2")),
        virtualSelectInput(paste0("error_type", sfx), "Error Type", choices = c("SD" = "sd", "SEM" = "sem"))
      ),
      layout_column_wrap(
        width = 1/2,
        checkboxInput(paste0("ShowErrorBar", sfx), "Error Bars", FALSE),
        checkboxInput(paste0("show_ref_line", sfx), "Ref. Line", FALSE)
      ),
      layout_column_wrap(
        width = 1/2,
        checkboxInput(paste0("show_significance", sfx), "Significance (*)", FALSE),
        checkboxInput(paste0("facet_by_gene", sfx), "Facet by Gene", FALSE)
      )
    ),
    accordion_panel(
      "Plot Appearance", icon = bsicons::bs_icon("palette"),
      textInput(paste0("plotTitle", sfx), "Plot Title"),
      layout_column_wrap(
        width = 1/2,
        textInput(paste0("xlab_custom", sfx), "X-Axis Label"),
        textInput(paste0("ylab_custom", sfx), "Y-Axis Label")
      ),
      virtualSelectInput(paste0("colorpicker", sfx), "Color Palette", choices = colorlist),
      layout_column_wrap(
        width = 1/2,
        sliderInput(paste0("axis_text_size", sfx), "Font Size", min = 6, max = 20, value = 12),
        sliderInput(paste0("xTextAngle", sfx), "X-Angle", min = 0, max = 90, value = 0)
      ),
      layout_column_wrap(
        width = 1/2,
        textInput(paste0("legendTitle", sfx), "Legend Title"),
        virtualSelectInput(paste0("legendPosition", sfx), "Position", choices = c("right", "none"))
      ),
      orderInput(paste0("x_order", sfx), "Sample Ordering", items = NULL, width = "100%"),
      span(class = "help-block", "Drag and drop to reorder.")
    ),
    accordion_panel(
      "Export", icon = bsicons::bs_icon("download"),
      layout_column_wrap(
        width = 1/2,
        numericInput(paste0("plot_w_", exp), "Width (in)", value = 7),
        numericInput(paste0("plot_h_", exp), "Height (in)", value = 5)
      ),
      layout_column_wrap(
        width = 1/2,
        virtualSelectInput(paste0("plot_format_", exp), "Format", choices = format_choices),
        numericInput(paste0("plot_res_", exp), "Scale", value = 3)
      ),
      downloadButton(download_id, "Download Plot", class = "btn-primary w-100")
    )
  )
}

####
# User Interface
####
ui <- page_navbar(
  title = "qRAT",
  id = "tabs",
  navbar_options = navbar_options(
    collapsible = TRUE),
  theme = bs_theme(
    bootswatch = "flatly",
    version = 5,
    primary = primary_col, 
    success = success_col,
    warning = warning_col,
    info = info_col,
    # Inter is served from www/fonts/inter (see header); no download from Google at startup
    base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif"),
    "navbar-bg" = primary_col,
    "body-bg" = page_bg_col
  ),
  
  header = tagList(
    # self-hosted Inter font (SIL Open Font License, see www/fonts/inter/OFL.txt)
    tags$link(rel = "stylesheet", href = "www/fonts/inter/font.css"),
    useWaiter(),
    waiterShowOnLoad(
      html = tagList(
        tags$div(
          style = sprintf("color: %s;", primary_col), 
          spin_fading_circles()
        ),
        tags$br(),
        tags$h4("Loading qRAT...", style = sprintf("color: %s;", primary_col))
      ),
      color = "#ffffff"
    ),
    
    
    shinyjs::useShinyjs(),
    
    
    tags$style(HTML(paste0("
      /* General Layout */
      body { background-color: ", page_bg_col, "; }
      .help-block { color: ", accent_blue, " !important; font-size: 0.85rem; margin-top: 5px; margin-bottom: 15px; }
      
      /* Card Design & Interactive Feedback */
      .card { 
        margin-bottom: 24px; 
        border-radius: 12px; 
        box-shadow: 0 4px 12px rgba(0,0,0,0.06); 
        border: 1px solid #dee2e6; 
        background-color: ", card_bg_col, " !important;
        transition: box-shadow 0.2s ease;
      }

      /* auto height only outside full-screen mode (full screen needs the card to fill the viewport) */
      .card:not([data-full-screen='true']) { height: auto !important; }

      /* no translate on hover: it shifted plots under the cursor and breaks full-screen positioning */
      .card:hover {
        box-shadow: 0 6px 18px rgba(0,0,0,0.1);
      }
      
      .card-header { 
        font-weight: 600; 
        background-color: ", card_header_bg, " !important;
        border-bottom: 1px solid #dcdfe3 !important;
        padding: 0.85rem 1.25rem;
        color: #2c3e50;
        display: flex;
        align-items: center;
        gap: 10px;
      }
      
      .card-body {
        padding: 1.25rem;
      }
      .card:not([data-full-screen='true']) .card-body {
        height: auto !important;
        overflow: visible !important;
      }

      /* top navbar: Flatly colours the active/hovered link green, hard to read on the blue bar */
      .navbar .nav-link.active, .navbar .nav-link:hover, .navbar .nav-link:focus { color: #ffffff !important; }
      .navbar .nav-link.active { box-shadow: inset 0 -3px 0 #ffffff; }

      /* ggiraph SVGs keep their aspect ratio: the container follows the SVG height; in full screen limit it to the viewport */
      .girafe.html-widget { height: auto !important; }
      .bslib-card[data-full-screen='true'] .girafe svg { max-height: calc(100vh - 160px); }

      /* main tab bar (Raw Data, Filtering, ...) */
      .nav-pills { margin-bottom: 1rem; row-gap: 0.25rem; }

      /* accordion in the sidebar: compact headers */
      .sidebar .accordion-button { font-weight: 600; font-size: 0.9rem; padding: 0.6rem 0.9rem; }
      
      /* spacing below the interactive plots */
      .girafe.html-widget {
        margin-bottom: 10px;
      }

      /* Sub-Tabs Styling */
      .navset-card-tab > .card-header {
        display: none;
      }
      
      /* Sidebar Elements */
      .sidebar h5 { font-weight: 700; margin-top: 20px; color: ", accent_blue, "; }
      .sidebar hr { margin: 15px 0; opacity: 0.1; }
      
      /* Button Styling */
      .btn-primary { 
        background-color: ", primary_col, "; 
        border-color: ", primary_col, "; 
        border-radius: 8px;
        font-weight: 500;
      }
      
      /* Navigation Tabs & Pills */
      .nav-pills .nav-link {
        border-radius: 8px;
        color: ", primary_col, ";
        font-weight: 500;
        margin-right: 5px;
      }
      
      .nav-pills .nav-link.active {
        background-color: ", accent_blue, " !important;
        color: white !important;
        box-shadow: 0 2px 6px rgba(50, 93, 136, 0.3);
      }
      
      .jumbotron-card { 
        background-color: #f8f9fa; 
        border-radius: 30px; 
        padding: 80px 40px; 
        border: 1px solid #e9ecef;
        margin-top: 20px; 
      }
    "))),
  
  
  # Einbinden des CSS-Stylings für die Sidebar
  tags$head(
    tags$style(HTML("
      .sidebar-section {
        background-color: #fcfcfc;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 16px;
        margin-bottom: 20px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.02);
      }
      .sidebar-section h6 {
        margin-top: 0;
        margin-bottom: 12px;
        font-size: 0.85rem;
        font-weight: 700;
        text-transform: uppercase;
        color: #475569;
        border-bottom: 1.5px solid #f1f5f9;
        padding-bottom: 8px;
      }
      .help-block {
        font-size: 0.75rem;
        color: #94a3b8;
        margin-top: 4px;
        line-height: 1.2;
        display: block;
      }
      .well { background-color: #ffffff !important; border: 1px solid #e2e8f0 !important; }
    "))
    )
  ),
  

  # --- HOME PAGE ---
  nav_panel(
    title = "Home",
    icon = icon("house"),
    div(class = "container py-5",
        # Jumbotron / Hero Section
        div(class = "jumbotron-card text-center mb-5",
            h1("qRAT", class = "display-3 fw-bold mb-3", style = paste0("color: ", accent_blue, ";")),
            p(class = "lead text-secondary mb-4", "qPCR Relative Expression Analysis Tool"),
            hr(class = "my-4", style = paste0("max-width: 100px; margin: auto; border-top: 4px solid ", primary_col, "; opacity: 1;")),
            p(class = "mx-auto mb-5", style = "max-width: 650px; font-size: 1.15rem; color: #666; line-height: 1.6;",
              "Analyze your qPCR data with precision and ease. qRAT provides an intuitive interface for both single and multi-plate workflows."
            ),
            div(class = "d-flex justify-content-center gap-3",
                actionButton("link_to_tab_single", "Single Plate Workflow", class = "btn btn-primary btn-lg px-5"),
                actionButton("link_to_tab_multi", "Multi Plate Workflow", class = "btn btn-outline-primary btn-lg px-5")
            )
        ),
        
        # Citation Section
        div(class = "row justify-content-center",
            div(class = "col-lg-10",
                div(class = "card border-0 shadow-sm", style = "background-color: #fcfcfc; border-radius: 15px;",
                    div(class = "card-body p-4",
                        div(class = "d-flex align-items-center mb-3",
                            icon("quote-left", class = "me-3", style = paste0("color: ", primary_col, "; font-size: 1.5rem;")),
                            h4("How to Cite", class = "mb-0", style = "font-weight: 700; color: #34495e;")
                        ),
                        p("If you use qRAT in published research, please support the developers by citing the following publications:", 
                          class = "text-muted mb-4"),
                        
                        div(class = "row g-4",
                            # Main Citation
                            div(class = "col-md-4",
                                div(style = "border-left: 3px solid #325d88; padding-left: 15px;",
                                    strong("for qRAT", style = "display: block; margin-bottom: 5px;"),
                                    tags$a(href = "https://doi.org/10.1186/s12859-022-04823-7", 
                                           "doi.org/10.1186/s12859-022-04823-7", 
                                           target = "_blank",
                                           style = "text-decoration: none; color: #325d88; font-size: 0.9rem;")
                                )
                            ),
                            # ddCt Citation
                            div(class = "col-md-4",
                                div(style = "border-left: 3px solid #5a7d9a; padding-left: 15px;",
                                    strong("for ddCq", style = "display: block; margin-bottom: 5px;"),
                                    tags$a(href = "https://doi.org/doi:10.18129/B9.bioc.ddCt", 
                                           "Bioconductor: ddCt", 
                                           target = "_blank",
                                           style = "text-decoration: none; color: #5a7d9a; font-size: 0.9rem;")
                                )
                            ),
                            # Limma Citation
                            div(class = "col-md-4",
                                div(style = "border-left: 3px solid #5a7d9a; padding-left: 15px;",
                                    strong("for statistics", style = "display: block; margin-bottom: 5px;"),
                                    tags$a(href = "https://doi.org/doi:10.18129/B9.bioc.limma", 
                                           "Bioconductor: limma", 
                                           target = "_blank",
                                           style = "text-decoration: none; color: #5a7d9a; font-size: 0.9rem;")
                                )
                            )
                        )
                    )
                )
            )
        )
        )
    ),
  
  
  # --- SINGLE PLATE ANALYSIS ---
  nav_panel(
    title = "Single Plate Analysis",
    icon = icon("stop"),
    layout_sidebar(
      fillable = FALSE, border = FALSE, border_radius = FALSE,
      sidebar = sidebar(
        width = 350,
        title = "Analysis Parameters",
        
        # --- 1. RAW DATA ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Raw Data'", 
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Import"), 
            radioGroupButtons(
              "Sep", "Data separator", 
              choices = c("auto", "comma", "tab", "semicolon"),
              status = "primary", size = "sm"
            ), 
            span(class="help-block", "Choose the separator used in the file. If unsure, use 'auto'."),
            fileInput("dtfile", "Upload (.csv, .txt)", width = "100%"),
            span(class="help-block", "Upload your raw Cq data file.")
          )
        ),
        
        # --- 2. FILTERING & QUALITY ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Filtering & Quality'", 
          # Experimental Controls
          tags$div(
            class = "sidebar-section",
            tags$h6("Experimental Controls"),
            virtualSelectInput("NTC_Input", "NTC/RT- Samples", choices = NULL, multiple = TRUE),
            span(class="help-block", "Exclude NTCs or RT- from calculations.")
          ),
          # --- FILTERING SECTION ---
          tags$div(
            class = "sidebar-section",
            tags$h6("Filtering Strategy"),
            
            # Auswahl des Modus
            radioButtons("filterMode", NULL, 
                         choices = c("Automatic Filtering" = "auto", 
                                     "Manual Filtering" = "manual"),
                         selected = "auto", inline = TRUE),
            br(),
            
            # --- AUTOMATIC FILTERING UI ---
            conditionalPanel(
              condition = "input.filterMode == 'auto'",
              tags$label("Cq Cut-Off", style="font-weight:500; font-size: 0.9rem;"),
              sliderInput("maxCt", NULL, min = 0, max = 45, value = c(5, 35)), 
              span(class="help-block", "Cq values outside this range will be excluded."),
              br(),
              tags$label("Replicate Variability", style="font-weight:500; font-size: 0.9rem;"),
              sliderInput("repDev", NULL, min = 0, max = 2, value = 0.3, step = 0.01),
              span(class="help-block", "Maximum allowed SD of the technical replicates. If exceeded, the replicate furthest from the mean is removed (at least 2 are kept).")
            ),
            
            # --- MANUAL FILTERING UI ---
            conditionalPanel(
              condition = "input.filterMode == 'manual'",
              p("Select specific wells/samples to exclude manually."),
              actionButton("open_manual_filter", "Select Outliers", icon = icon("list"), class = "btn-primary btn-sm w-100"),
              br(), br(),
              uiOutput("manual_filter_status"), # Zeigt an, wie viele ausgeschlossen wurden
              span(class="help-block", "Data points marked in the list will be removed from analysis.")
            )
          )
        ),
        
        
        
        conditionalPanel(
          condition = "input.tabsSingle=='Reference Finder'",
          tags$div(
            class = "sidebar-section",
            tags$h6("Stability Analysis"),
              # Auswahl der Gene, die getestet werden sollen
              virtualSelectInput(
                "ref_candidates", "Select Candidate Genes",
                choices = NULL, # Wird dynamisch gefüllt
                multiple = TRUE,
                search = TRUE,
                placeholder = "Select at least 3 genes..."
              ),
              hr(),
              tags$p(class = "small text-muted mb-2", "Algorithm: geNorm (M value)"),
              actionButton("run_val_btn", "Run Validation", 
                           class = "btn-primary w-100", icon = icon("play")),
              hr(),
              # Quick-Action
              helpText("Once validated, you can apply the best genes to your analysis."),
              actionButton("apply_best_ref", "Apply Top Genes as Refs", 
                           class = "btn-outline-success btn-sm w-100")
          )
        ),
            
            
        
        
        
        
        # ---  RELATIVE dCq and ddCq ---
        conditionalPanel(
          condition = "input.tabsSingle=='Relative dCq' || input.tabsSingle=='Relative ddCq'",
          # Data Selection
          tags$div(
            class = "sidebar-section",
            tags$h6("Normalization & Grouping"),
            virtualSelectInput("Refs", "Reference Gene(s)", choices = "", multiple = TRUE, placeholder = "Select..."),
            hr(),
            checkboxInput("groupSamples", "Group Biological Replicates", FALSE),
            conditionalPanel(
              condition = "input.groupSamples == true",
              actionButton("open_grouping_modal", "Configure Groups", class = "btn-primary btn-sm w-100"),
              p(class="help-block", "Define groups and assign samples. Each sample then counts as one biological replicate (mean of its technical replicates).")
            ),
            hr(),
        )),
        

        # --- 3. RELATIVE dCq ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Relative dCq'",
          plot_settings_sidebar("", "dCq", "download_plot_dCq", c("RQ", "dCq", "-dCq"), c("pdf", "svg", "png"))
        ),
        
        
        # --- 4. RELATIVE ddCq ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Relative ddCq'",
          plot_settings_sidebar("DDCt", "ddCq", "download_plot_ddCq", c("Fold Change", "ddCt"), c("png", "svg", "pdf"), calibrator = tagList(virtualSelectInput("Mock", "Calibrator (Control)", choices = "", placeholder = "Select baseline..."), span(class="help-block", "The baseline sample for ddCq."), hr()))
        ),
        
        # --- 5. STATISTICAL ANALYSIS ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Statistics'",
          # Analysis Input
          tags$div(
            class = "sidebar-section",
            tags$h6("Analysis Input"),
            virtualSelectInput("limma_input", "Linear Model Input", choices = c("dCq"))
          ),
          # Parameters
          tags$div(
            class = "sidebar-section",
            tags$h6("Parameters"),
            virtualSelectInput("adjustMethod", "P-Value Adjustment", choices = c("BH", "Holm", "Bonferroni")),
            br(),
            radioButtons("compType", "Comparison Type", choices = c("Single Comparison", "Multiple paired Comparisons")),
            hr(),
            # Comparison Setup
            conditionalPanel(
              condition = "input.compType == 'Single Comparison'",
              virtualSelectInput("Comps", "Comparison Reference", choices = "", placeholder = "Choose control..."),
              span(class="help-block", "Compare all samples against this control.")
            ),
            conditionalPanel(
              condition = "input.compType == 'Multiple paired Comparisons'",
              tags$label("Paired Comparisons", style="font-weight:bold;"),
              layout_column_wrap(
                width = 1/2,
                actionButton("add_btn", "Add", class = "btn-primary btn-sm w-100"),
                actionButton("rm_btn", "Remove", class = "btn-outline-danger btn-sm w-100")
              ),
              uiOutput("comparisonBoxes_ui")
            )
          )
        )
      ), # Ende Sidebar
      
      
        # Hauptbereich Single Plate Analyse mit Tabs
        navset_pill(
          id = "tabsSingle",
          nav_panel("Raw Data",
                    icon = icon("table"),
                    layout_column_wrap(
                      width = 1/2,
                      card(full_screen = TRUE, card_header(tagList(icon("table"), "Imported Data Overview")), card_body(DTOutput("dataSinglePlate") %>% withSpinner())),
                      navset_card_underline(
                        full_screen = TRUE,
                        nav_panel("Plate View", icon = icon("eye"), girafeOutput("plotCtCard", height = "500px") %>% withSpinner(), plate_selection_bar("apply_plate_selection")),
                        nav_panel("Global Distribution", icon = icon("chart-area"), girafeOutput("plotCtDistrib", height = "500px") %>% withSpinner())
                      )
                    )
          ),
          nav_panel("Filtering & Quality",
                    icon = icon("vial-circle-check"),
                    layout_column_wrap(
                      width = 1/2,
                      card(full_screen = TRUE, card_header(tagList(icon("filter"), "Well Exclusion Review")),
                           card_body(
                             p(class="text-muted small mb-2", "Wells excluded by the current filter settings. Switch to 'Manual Filtering' in the sidebar to choose exclusions yourself."),
                             DTOutput("dataSinglePlateBadRep") %>% withSpinner()
                           )),
                      navset_card_underline(
                        full_screen = TRUE,
                        nav_panel("Outlier Analysis", girafeOutput("SinglePlateFilterPlot", height = "500px")),
                        nav_panel("Cq Distribution", girafeOutput("SinglePlateBoxplot", height = "500px")),
                        nav_panel("MIQE Compliance", uiOutput("SinglePlateMIQEcheck"))
                      )
                    )
          ),
          nav_panel("Reference Finder",
                    icon = icon("list-check"),
                    layout_column_wrap(
                      width = 1/2,
                      navset_card_underline(
                        full_screen = TRUE,
                        nav_panel("Ranking Table", DTOutput("ref_ranking_table")),
                        nav_panel("Stability Plot", girafeOutput("ref_stability_plot"))
                      ),
                      card(
                        fill = FALSE,
                        card_header("Algorithm Info"),
                        card_body(
                          uiOutput("algo_description")
                        )
                      ))
          ),

          nav_panel("Relative dCq",
                    icon = icon("chart-bar"),
                    card(fill = FALSE, full_screen = TRUE, card_header(tagList(icon("chart-bar"), "Gene Expression (dCq)")), card_body(girafeOutput("dcqPlot", height = "500px") %>% withSpinner())),
                    card(fill = FALSE, full_screen = TRUE, card_header(tagList(icon("list"), "Normalized Data Table")), card_body(DTOutput("ddctAbsolute")))
          ),
          nav_panel("Relative ddCq",
                    icon = icon("chart-line"),
                    card(fill = FALSE, full_screen = TRUE, card_header(tagList(icon("chart-line"), "Relative Expression (Fold Change)")), card_body(girafeOutput("ddctRelGraph2", height = "500px") %>% withSpinner())),
                    card(fill = FALSE, full_screen = TRUE, card_header(tagList(icon("clipboard-list"), "ddCq Details")), card_body(DTOutput("ddctRelative")))
          ),


          nav_panel("Statistics",
                    icon = icon("calculator"),
                    card(full_screen = TRUE, card_header(tagList(icon("calculator"), "Limma Statistics Results")), card_body(DTOutput("resultLimma") %>% withSpinner()))
          )
        )

        )# Ende layout_sidebar
      ), # Ende Single Plate Analysis nav_panel
  
  
  # --- MULTIPLE PLATE ANALYSIS ---
  nav_panel(
    title = "Multiple Plate Analysis",
    icon = icon("table-cells-large"),
    layout_sidebar(
      fillable = FALSE, border = FALSE, border_radius = FALSE,
      sidebar = sidebar(
        width = 350,
        title = "Analysis Parameters",
        
        # --- 1. RAW DATA MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Raw Data'",
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Import"),
            radioGroupButtons(
              "SepM", "Data separator", 
              choices = c("auto", "comma", "tab", "semicolon"), 
              status = "primary", size = "sm"
            ),
            span(class="help-block", "Choose the separator used in the files."),
            fileInput("plates", "Upload multiple plates (.csv, .txt)", multiple = TRUE, width = "100%"),
            span(class="help-block", "Hold [Ctrl] or [Cmd] to select multiple files.")
          )
        ),
        
        # --- 2. FILTERING & QUALITY MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Filtering & Quality'",
          # Experimental Controls
          tags$div(
            class = "sidebar-section",
            tags$h6("Experimental Controls"),
            virtualSelectInput("NTC_Input_MP", "NTC/RT- Samples", choices = NULL, multiple = TRUE),
            span(class="help-block", "Exclude NTCs or RT- from calculations.")
          ),
          # --- FILTERING STRATEGY MULTI PLATE ---
          tags$div(
            class = "sidebar-section",
            tags$h6("Filtering Strategy"),
            
            # Auswahl des Modus
            radioButtons("filterModeMP", NULL, 
                         choices = c("Automatic Filtering" = "auto", 
                                     "Manual Filtering" = "manual"),
                         selected = "auto", inline = TRUE),
            br(),
            
            # --- AUTOMATIC FILTERING UI (MP) ---
            conditionalPanel(
              condition = "input.filterModeMP == 'auto'",
              tags$label("Cq Cut-Off", style="font-weight:500; font-size: 0.9rem;"),
              sliderInput("maxCtMulti", NULL, min = 0, max = 45, value = c(5, 35)),
              span(class="help-block", "Cq values outside this range will be excluded."),
              br(),
              tags$label("Replicate Variability", style="font-weight:500; font-size: 0.9rem;"),
              sliderInput("repDevMulti", NULL, min = 0, max = 2, value = 0.3, step = 0.01),
              span(class="help-block", "Maximum allowed SD of the technical replicates. If exceeded, the replicate furthest from the mean is removed (at least 2 are kept).")
            ),
            
            # --- MANUAL FILTERING UI (MP) ---
            conditionalPanel(
              condition = "input.filterModeMP == 'manual'",
              p("Select specific wells/samples to exclude manually."),
              actionButton("open_manual_filter_MP", "Select Outliers", icon = icon("list"), class = "btn-primary btn-sm w-100"),
              br(), br(),
              uiOutput("manual_filter_status_MP"), # Status Anzeige für MP
              span(class="help-block", "Data points marked in the list will be removed from analysis.")
            )
          )
        ),
        
        # --- 3. INTER-PLATE CALIBRATION ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Inter-Plate Calibration'",
          tags$div(
            class = "sidebar-section",
            tags$h6("Calibration Calculation"),
            prettySwitch("Id027", "Use Calibration", status = "success", fill = TRUE, slim = FALSE),
            virtualSelectInput("IPC", "Select IPC Sample", choices = "", placeholder = "Select..."),
            span(class="help-block", "Choose the sample present on all plates.")
          ),
          tags$div(
            class = "sidebar-section",
            tags$h6("Calibration Check"),
            virtualSelectInput("SamplePickerIPCcomparison", "Samples to Plot", choices = "", multiple = TRUE, search = TRUE),
            span(class="help-block", "Visualize Cq stability across plates.")
          )
        ),
        
        
        
        # --- REFERENCE FINDER SIDEBAR MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Reference Finder'",
          tags$div(
            class = "sidebar-section",
            tags$h6("Stability Analysis (Multi-Plate)"),
            # Auswahl der Kandidaten-Gene für Multi-Plate
            virtualSelectInput(
              "ref_candidates_multi", "Select Candidate Genes",
              choices = NULL, # Wird dynamisch über multiData2() befüllt
              multiple = TRUE,
              search = TRUE,
              placeholder = "Select at least 3 genes..."
            ),
            hr(),
            tags$p(class = "small text-muted mb-2", "Algorithm: geNorm (M value)"),
            actionButton("run_val_btn_multi", "Run Validation", 
                         class = "btn-primary w-100", icon = icon("play")),
            hr(),
            # Quick-Action für Multi-Plate
            helpText("Apply the most stable genes to your Multi-Plate analysis."),
            actionButton("apply_best_ref_multi", "Apply Top Genes as Refs", 
                         class = "btn-outline-success btn-sm w-100")
          )
        ),
        
        
        # ---  RELATIVE dCq and ddCq Multi---
        conditionalPanel(
          condition = "input.tabsMulti=='Relative dCq' || input.tabsMulti=='Relative ddCq'",
          # Data Selection
          tags$div(
            class = "sidebar-section",
            tags$h6("Normalization & Grouping"),
            virtualSelectInput("RefsM", "Reference Gene(s)", choices = "", multiple = TRUE, placeholder = "Select..."),
            hr(),
            checkboxInput("groupSamplesMulti", "Group Biological Replicates", FALSE),
            conditionalPanel(
              condition = "input.groupSamplesMulti == true",
              actionButton("open_grouping_modalMulti", "Configure Groups", class = "btn-primary btn-sm w-100"),
              p(class="help-block", "Define groups and assign samples. Each sample then counts as one biological replicate (mean of its technical replicates).")
            ),
            hr(),
          )),
        
        
        # --- 4. RELATIVE dCq MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Relative dCq'",
          plot_settings_sidebar("Multi", "Multi", "download_plot_MultiAbs", c("RQ", "dCq", "-dCq"), c("pdf", "svg", "png"))
        ),
        
        # --- 5. RELATIVE ddCq MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Relative ddCq'",
          plot_settings_sidebar("DDCtMulti", "DDCtMulti", "download_plot_MultiRel", c("Fold Change", "ddCt"), c("png", "svg", "pdf"), calibrator = tagList(virtualSelectInput("MockM", "Calibrator (Control)", choices = NULL, placeholder = "Select baseline..."), span(class="help-block", "The baseline sample for ddCq."), hr()))
        ),
        
        # --- 6. STATISTICAL ANALYSIS MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Statistics'",
          # Analysis Input
          tags$div(
            class = "sidebar-section",
            tags$h6("Analysis Input"),
            virtualSelectInput("limma_inputMulti", "Linear Model Input", choices = c("dCq"))
          ),
          # Parameters
          tags$div(
            class = "sidebar-section",
            tags$h6("Parameters"),
            virtualSelectInput("adjustMethodMulti", "P-Value Adjustment", choices = c("BH", "Holm", "Bonferroni")),
            br(),
            radioButtons("compTypeM", "Comparison Type", choices = c("Single Comparison", "Multiple paired Comparisons")),
            hr(),
            # Comparison Setup
            conditionalPanel(
              condition = "input.compTypeM == 'Single Comparison'",
              virtualSelectInput("CompsM", "Comparison Reference", choices = "", placeholder = "Choose control..."),
              span(class="help-block", "Compare all samples against this control.")
            ),
            conditionalPanel(
              condition = "input.compTypeM == 'Multiple paired Comparisons'",
              tags$label("Paired Comparisons", style="font-weight:bold;"),
              layout_column_wrap(
                width = 1/2,
                actionButton("add_btnM", "Add", class = "btn-primary btn-sm w-100"),
                actionButton("rm_btnM", "Remove", class = "btn-outline-danger btn-sm w-100")
              ),
              uiOutput("comparisonBoxesM_ui")
            )
          )
        )
      ),
      
      # Hauptbereich Multi Plate
      navset_pill(
        id = "tabsMulti",

        # --- TAB 1: RAW DATA ---
        nav_panel(
          title = "Raw Data",
          icon = icon("table"),
          conditionalPanel(
            condition = "output.fileUploadedMulti",
            layout_column_wrap(
              width = 1/2,
              card(
                full_screen = TRUE,
                card_header(tagList(icon("table"), "Imported Multi-Plate Data")),
                card_body(DTOutput("multiplePlatesData") %>% withSpinner())
              ),
              navset_card_underline(
                full_screen = TRUE,
                nav_panel(
                  title = "Plate View",
                  icon = icon("eye"),
                  # Plattenauswahl als Toolbar über dem Plot
                  div(
                    class = "d-flex justify-content-between align-items-center gap-2",
                    tags$span("Cq per Plate", class = "fw-semibold"),
                    virtualSelectInput(
                      "PlateSelect",
                      label = NULL,
                      choices = "",
                      placeholder = "Select Plate",
                      width = "180px"
                    )
                  ),
                  girafeOutput("MultiplotCtCard", height = "500px"),
                  plate_selection_bar("apply_plate_selection_MP", "Selection applies to the plate shown; exclusions on other plates are kept.")
                ),
                nav_panel(
                  title = "Global Distribution",
                  icon = icon("chart-area"),
                  girafeOutput("MultiplotCtDistrib", height = "500px")
                )
              )
            )
          )
        ),

        # --- TAB 2: FILTERING & QUALITY ---
        nav_panel(
          title = "Filtering & Quality",
          icon = icon("vial-circle-check"),
          conditionalPanel(
            condition = "output.fileUploadedMulti",
            layout_column_wrap(
              width = 1/2,
              card(
                full_screen = TRUE,
                card_header(tagList(icon("filter"), "Well Exclusion Review")),
                card_body(
                  p(class="text-muted small mb-2", "Wells excluded by the current filter settings (all plates). Switch to 'Manual Filtering' in the sidebar to choose exclusions yourself."),
                  DTOutput("dataMultiplePlatesBadRep") %>% withSpinner()
                )
              ),
              navset_card_underline(
                full_screen = TRUE,
                nav_panel("Outlier Analysis", girafeOutput("MultiplePlatesFilterPlot", height = "500px")),
                nav_panel("Cq Distribution", girafeOutput("MultiplePlatesBoxplot", height = "500px")),
                nav_panel("MIQE Compliance", uiOutput("MultiplePlateMIQEcheck"))
              )
            )
          )
        ),

        # --- TAB 3: INTER-PLATE CALIBRATION ---
        nav_panel(
          title = "Inter-Plate Calibration",
          icon = icon("scale-balanced"),
          conditionalPanel(
            condition = "output.fileUploadedMulti",

            # 2:1 - links die Analyse, rechts die Methodik (auf kleinen Bildschirmen untereinander)
            layout_columns(
              col_widths = c(8, 4),
              navset_card_underline(
                full_screen = TRUE,
                # Tab 1: Die Roh-IPCs
                nav_panel(
                  title = "Inter-Plate Calibrators",
                  DTOutput("extractedIPCsTable")
                ),
                # Tab 2: Der visuelle Check
                nav_panel(
                  title = "Calibration Check Plot",
                  girafeOutput("plotIPCcomparison", height = "550px")
                ),
                # Tab 3: Die berechneten Faktoren
                nav_panel(
                  title = "Calibration Factors",
                  DTOutput("tableCalibrationFactors") %>% withSpinner()
                )
              ),
              card(
                card_header(tagList(icon("info-circle"), "Methodology")),
                card_body(
                  uiOutput("ipc_info_text")
                )
              )
            ),

            # UNTERER BEREICH: Die resultierenden Daten
            card(
              full_screen = TRUE,
              card_header(tagList(icon("table"), "Full Calibrated Data Table")),
              card_body(
                DTOutput("interPlateCalibration")
              )
            )
          )
        ),


        # --- REFERENCE FINDER MAIN MULTI ---
        nav_panel("Reference Finder",
                  icon = icon("list-check"),
                  layout_column_wrap(
                    width = 1/2,
                    # Eigene IDs für die Multi-Plate Outputs
                    navset_card_underline(
                      full_screen = TRUE,
                      nav_panel("Ranking Table", DTOutput("ref_ranking_table_multi")),
                      nav_panel("Stability Plot", girafeOutput("ref_stability_plot_multi"))
                    ),
                    card(
                      fill = FALSE,
                      card_header("Algorithm Info"),
                      card_body(
                        uiOutput("algo_description_multi")
                      )
                    )
                  )
        ),


        # --- TAB 4: RELATIVE dCq ---
        nav_panel(
          title = "Relative dCq",
          icon = icon("chart-column"),
            card(fill = FALSE, full_screen = TRUE,
              card_header(tagList(icon("chart-bar"), "Gene Expression (dCq)")),
              card_body(girafeOutput("ddctMultiAbsGraph", height = "500px") %>% withSpinner())
            ),
            card(fill = FALSE, full_screen = TRUE,
              card_header(tagList(icon("list"), "Normalized Data Table")),
              card_body(DTOutput("ddctAbsoluteMulti"))
            )
        ),

        # --- TAB 5: RELATIVE ddCq ---
        nav_panel(
          title = "Relative ddCq",
          icon = icon("chart-line"),
            card(fill = FALSE, full_screen = TRUE,
              card_header(tagList(icon("chart-line"), "Relative Expression (Fold Change)")),
              card_body(girafeOutput("ddctRelGraphMulti", height = "500px") %>% withSpinner())
            ),
            card(fill = FALSE, full_screen = TRUE,
              card_header(tagList(icon("clipboard-list"), "ddCq Details")),
              card_body(DTOutput("ddctRelativeMulti"))
            )
        ),

        # --- TAB 6: STATISTICAL ANALYSIS ---
        nav_panel(
          title = "Statistics",
          icon = icon("calculator"),
          card(
            full_screen = TRUE,
            card_header(tagList(icon("calculator"), "Limma Statistics Results")),
            card_body(
              p("Statistical analysis of differentially expressed genes for multi-plate experiments."),
              DTOutput("resultLimmaMulti") %>% withSpinner()
            )
          )
        )
      )
    )
  ),
    
    # --- GLOBAL HELP (In the Header) ---
    nav_panel(
      title = "Help & About", 
      icon = icon("circle-question"), 
      includeHTML("help.html")
    )
    ) # Ende page_navbar
