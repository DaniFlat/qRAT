#' @name qRAT - qPCR Relative Expression Analysis Tool
#' @title UI for the qRAT Shiny app
#' @author Daniel Flatschacher

library("shiny")
library("bslib")
library("shinyWidgets")
library("plotly")
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
    base_font = font_google("Inter"),
    "navbar-bg" = primary_col,
    "body-bg" = page_bg_col
  ),
  
  header = tagList(
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
        height: auto !important; 
        transition: transform 0.2s ease, box-shadow 0.2s ease;
      }
      
      .card:hover {
        transform: translateY(-2px);
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
        height: auto !important; 
        overflow: visible !important; 
      }
      
      /* Ensure Plotly containers take up required space without scrolling */
      .js-plotly-plot, .plotly {
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
    ")))
  ),
  
  
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
    page_sidebar(
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
              span(class="help-block", "Maximum allowed standard deviation between replicates.")
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
              # Parameter für die Algorithmen
              tags$h6("Parameters"),
              selectInput("norm_method", "Primary Algorithm", 
                          choices = c("NormFinder", "geNorm")),
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
              p(class="help-block", "Define groups and assign samples.")
            ),
            hr(),
        )),
        

        # --- 3. RELATIVE dCq ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Relative dCq'",
          # Data Selection
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Selection"),
            hr(),
            virtualSelectInput("SamplePicker", "Filter Samples", choices = NULL, multiple = TRUE, search = TRUE),
            virtualSelectInput("GenePicker", "Filter Genes", choices = NULL, multiple = TRUE, search = TRUE)
          ),
          # Data Visualization
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Visualization"),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("PlotData", "Data Type", choices = c("RQ", "dCq", "-dCq")),
              virtualSelectInput("PlotType", "Plot Type", choices = c("Bar Chart", "Dot Plot"))
            ),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("scale", "Y-Axis", choices = c("normal", "log10", "log2")),
              virtualSelectInput("error_type", "Error Type", choices = c("SD" = "sd", "SEM" = "sem"))
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("ShowErrorBar", "Error Bars", FALSE),
              checkboxInput("show_ref_line", "Ref. Line", FALSE)
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("show_significance", "Significance (*)", FALSE),
              checkboxInput("facet_by_gene", "Facet by Gene", FALSE)
            )
          ),
          # Plot Appearance
          switchInput("show_plot_settings_dCq", label = "Appearance", value = FALSE, size = "small", onStatus = "success"),
          conditionalPanel(
            condition = "input.show_plot_settings_dCq",
            tags$div(
              class = "sidebar-section",
              tags$h6("Plot Appearance"),
              textInput("plotTitle", "Plot Title"),
              layout_column_wrap(
                width = 1/2,
                textInput("xlab_custom", "X-Axis Label"),
                textInput("ylab_custom", "Y-Axis Label")
              ),
              virtualSelectInput("colorpicker", "Color Palette", choices = colorlist),
              layout_column_wrap(
                width = 1/2,
                sliderInput("axis_text_size", "Font Size", min = 6, max = 20, value = 12),
                sliderInput("xTextAngle", "X-Angle", min = 0, max = 90, value = 0)
              ),
              layout_column_wrap(
                width = 1/2,
                textInput("legendTitle", "Legend Title"),
                virtualSelectInput("legendPosition", "Position", choices = c("right", "none"))
              ),
              tags$label("Sample Ordering", style="font-size:0.8rem; font-weight:bold;"),
              orderInput("x_order", NULL, items = NULL, width = "100%"),
              span(class="help-block", "Drag and drop to reorder.")
            )
          ),
          # Export Options
          switchInput("show_export_settings_dCq", label = "Export Opt.", value = FALSE, size = "small"),
          conditionalPanel(
            condition = "input.show_export_settings_dCq",
            tags$div(
              class = "sidebar-section",
              tags$h6("Export Options"),
              layout_column_wrap(
                width = 1/2,
                numericInput("plot_w_dCq", "Width (in)", value = 7),
                numericInput("plot_h_dCq", "Height (in)", value = 5)
              ),
              layout_column_wrap(
                width = 1/2,
                virtualSelectInput("plot_format_dCq", "Format", choices = c("pdf", "svg", "png")),
                numericInput("plot_res_dCq", "Scale", value = 3)
              ),
              downloadButton("download_plot_dCq", "Download Plot", class = "btn-primary w-100")
            )
          )
        ),
        
        
        # --- 4. RELATIVE ddCq ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Relative ddCq'",
          # Data Selection (inkl. Calibrator)
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Selection"),
            virtualSelectInput("Mock", "Calibrator (Control)", choices = "", placeholder = "Select baseline..."),
            span(class="help-block", "The baseline sample for ddCq."),
            hr(),
            virtualSelectInput("SamplePickerDDCt", "Filter Samples", choices = "", multiple = TRUE),
            virtualSelectInput("GenePickerDDCt", "Filter Genes", choices = "", multiple = TRUE)
          ),
          # Data Visualization
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Visualization"),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("PlotDataDDCt", "Data Type", choices = c("Fold Change", "ddCt")),
              virtualSelectInput("PlotTypeDDCt", "Plot Type", choices = c("Bar Chart", "Dot Plot"))
            ),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("scaleDDCt", "Y-Axis", choices = c("normal", "log10", "log2")),
              virtualSelectInput("error_typeDDCt", "Error Value", choices = c("SD" = "sd", "SEM" = "sem"))
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("ShowErrorBarDDCt", "Error Bars", FALSE),
              checkboxInput("show_ref_lineDDCt", "Ref. Line", FALSE)
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("show_significanceDDCt", "Significance (*)", FALSE),
              checkboxInput("facet_by_geneDDCt", "Facet by Gene", FALSE)
            )
          ),
          # Appearance
          switchInput("show_appearance_settings_ddCq", label = "Appearance", value = FALSE, size = "small", onStatus = "success"),
          conditionalPanel(
            condition = "input.show_appearance_settings_ddCq",
            tags$div(
              class = "sidebar-section",
              tags$h6("Plot Appearance"),
              textInput("plotTitleDDCt", "Plot Title"),
              layout_column_wrap(
                width = 1/2,
                textInput("xlab_customDDCt", "X-Axis Label"),
                textInput("ylab_customDDCt", "Y-Axis Label")
              ),
              virtualSelectInput("colorpickerDDCt", "Color Palette", choices = colorlist),
              layout_column_wrap(
                width = 1/2,
                sliderInput("axis_text_sizeDDCt", "Font Size", min = 6, max = 20, value = 12),
                sliderInput("xTextAngleDDCt", "X-Angle", min = 0, max = 90, value = 0)
              ),
              layout_column_wrap(
                width = 1/2,
                textInput("legendTitleDDCt", "Legend Title"),
                virtualSelectInput("legendPositionDDCt", "Position", choices = c("right", "none"))
              ),
              orderInput("x_orderDDCt", "Sample Ordering", items = NULL, width = "100%"),
              span(class="help-block", "Drag and drop to reorder.")
            )
          ),
          # Export
          switchInput("show_export_settingsDDCt", label = "Export Opt.", value = FALSE, size = "small"),
          conditionalPanel(
            condition = "input.show_export_settingsDDCt",
            tags$div(
              class = "sidebar-section",
              tags$h6("Export Options"),
              layout_column_wrap(
                width = 1/2,
                numericInput("plot_w_ddCq", "Width (in)", value = 7), 
                numericInput("plot_h_ddCq", "Height (in)", value = 5)
              ),
              layout_column_wrap(
                width = 1/2,
                virtualSelectInput("plot_format_ddCq", "Format", choices = c("png", "svg", "pdf")),
                numericInput("plot_res_ddCq", "Scale", value = 3)
              ),
              downloadButton("download_plot_ddCq", "Download Plot", class = "btn-primary w-100")
            )
          )
        ),
        
        # --- 5. STATISTICAL ANALYSIS ---
        conditionalPanel(
          condition = "input.tabsSingle == 'Statistics'",
          # Analysis Input
          tags$div(
            class = "sidebar-section",
            tags$h6("Analysis Input"),
            virtualSelectInput("limma_input", "Linear Model Input", choices = c("dCq")),
            tags$div(
              style = "font-size: 0.8rem; padding: 10px; background: #fff3cd; color: #856404; border-radius: 4px; margin-top:10px;",
              tags$b("R Syntax:"), " Samples must start with a letter. Avoid special characters."
            )
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
        navset_card_pill(
          id = "tabsSingle",
          nav_panel("Raw Data",
                    icon = icon("table"),
                    layout_column_wrap(
                      width = 1/2,
                      card(card_header(tagList(icon("table"), "Imported Data Overview")), card_body(DTOutput("dataSinglePlate") %>% withSpinner())),
                      navset_card_underline(
                        nav_panel("Plate View", icon = icon("eye"), plotlyOutput("plotCtCard", height = "500px") %>% withSpinner()),
                        nav_panel("Global Distribution", icon = icon("chart-area"), plotlyOutput("plotCtDistrib", height = "500px") %>% withSpinner())
                      )
                    )
          ),
          nav_panel("Filtering & Quality",
                    icon = icon("vial-circle-check"),
                    layout_column_wrap(
                      width = 1/2,
                      card(card_header(tagList(icon("filter"), "Well Exclusion Review")), 
                           card_body(
                             p(class="text-muted small mb-2", "Review flagged wells. Uncheck to manually re-include a well or check to mark it as an outlier."),
                             DTOutput("dataSinglePlateBadRep") %>% withSpinner()
                           )),
                      navset_card_underline(
                        nav_panel("Outlier Analysis", plotlyOutput("SinglePlateFilterPlot", height = "500px")),
                        nav_panel("Cq Distribution", plotlyOutput("SinglePlateBoxplot", height = "500px")),
                        nav_panel("MIQE Compliance", uiOutput("SinglePlateMIQEcheck"))
                      )
                    )
          ),
          nav_panel("Reference Finder",
                    icon = icon("list-check"),
                    layout_column_wrap(
                      width = 1/2,
                          # Wir nutzen Tabs innerhalb der Card für verschiedene Ansichten
                          navset_card_underline(
                            nav_panel("Ranking Table", DTOutput("ref_ranking_table")),
                            nav_panel("Stability Plot", plotlyOutput("ref_stability_plot"))
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
                    card(fill = FALSE, card_header(tagList(icon("chart-bar"), "Gene Expression (dCq)")), card_body(plotlyOutput("dcqPlot", height = "500px") %>% withSpinner())),
                    card(fill = FALSE, card_header(tagList(icon("list"), "Normalized Data Table")), card_body(DTOutput("ddctAbsolute"))),
          ),
          nav_panel("Relative ddCq",
                    icon = icon("chart-line"),
                    card(fill = FALSE, card_header(tagList(icon("chart-line"), "Relative Expression (Fold Change)")), card_body(plotlyOutput("ddctRelGraph2", height = "500px") %>% withSpinner())),
                    card(fill = FALSE, card_header(tagList(icon("clipboard-list"), "ddCq Details")), card_body(DTOutput("ddctRelative")))
          ),
      
          
          nav_panel("Statistics",
                    icon = icon("calculator"),
                    card(card_header(tagList(icon("calculator"), "Limma Statistics Results")), card_body(DTOutput("resultLimma") %>% withSpinner()))
          )
        )

        )# Ende page_sidebar
      ), # Ende Single Plate Analysis nav_panel
  
  
  # --- MULTIPLE PLATE ANALYSIS ---
  nav_panel(
    title = "Multiple Plate Analysis",
    icon = icon("table-cells-large"),
    page_sidebar(
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
              span(class="help-block", "Maximum allowed standard deviation between replicates.")
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
            tags$h6("Parameters"),
            selectInput("norm_method_multi", "Primary Algorithm", 
                        choices = c("NormFinder")), # geNorm optional
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
              p(class="help-block", "Define groups and assign samples.")
            ),
            hr(),
          )),
        
        
        # --- 4. RELATIVE dCq MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Relative dCq'",
          # Data Selection
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Selection"),
          #  virtualSelectInput("RefsM", "Reference Gene(s)", choices = "", multiple = TRUE, placeholder = "Select..."),
            hr(),
            virtualSelectInput("SamplePickerMulti", "Filter Samples", choices = NULL, multiple = TRUE, search = TRUE),
            virtualSelectInput("GenePickerMulti", "Filter Genes", choices = NULL, multiple = TRUE, search = TRUE)
          ),
          # Data Visualization
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Visualization"),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("PlotDataMulti", "Data Type", choices = c("RQ", "dCq", "-dCq")),
              virtualSelectInput("PlotTypeMulti", "Plot Type", choices = c("Bar Chart", "Dot Plot"))
            ),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("scaleMulti", "Y-Axis", choices = c("normal", "log10", "log2")),
              virtualSelectInput("error_typeMulti", "Error Type", choices = c("SD" = "sd", "SEM" = "sem"))
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("ShowErrorBarMulti", "Error Bars", FALSE),
              checkboxInput("show_ref_lineMulti", "Ref. Line", FALSE)
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("show_significanceMulti", "Significance (*)", FALSE),
              checkboxInput("facet_by_geneMulti", "Facet by Gene", FALSE)
            )
          ),
          # Appearance
          switchInput("show_plot_settings_Multi", label = "Appearance", value = FALSE, size = "small", onStatus = "success"),
          conditionalPanel(
            condition = "input.show_plot_settings_Multi",
            tags$div(
              class = "sidebar-section",
              tags$h6("Plot Appearance"),
              textInput("plotTitleMulti", "Plot Title"),
              layout_column_wrap(
                width = 1/2,
                textInput("xlab_customMulti", "X-Axis Label"),
                textInput("ylab_customMulti", "Y-Axis Label")
              ),
              virtualSelectInput("colorpickerMulti", "Color Palette", choices = colorlist),
              layout_column_wrap(
                width = 1/2,
                sliderInput("axis_text_sizeMulti", "Font Size", min = 6, max = 20, value = 12),
                sliderInput("xTextAngleMulti", "X-Angle", min = 0, max = 90, value = 0)
              ),
              layout_column_wrap(
                width = 1/2,
                textInput("legendTitleMulti", "Legend Title"),
                virtualSelectInput("legendPositionMulti", "Position", choices = c("right", "none"))
              ),
              tags$label("Sample Ordering", style="font-size:0.8rem; font-weight:bold;"),
              orderInput("x_orderMulti", NULL, items = NULL, width = "100%"),
              span(class="help-block", "Drag and drop to reorder.")
            )
          ),
          # Export
          switchInput("show_export_settings_Multi", label = "Export Opt.", value = FALSE, size = "small"),
          conditionalPanel(
            condition = "input.show_export_settings_Multi",
            tags$div(
              class = "sidebar-section",
              tags$h6("Export Options"),
              layout_column_wrap(
                width = 1/2,
                numericInput("plot_w_Multi", "Width (in)", value = 7),
                numericInput("plot_h_Multi", "Height (in)", value = 5)
              ),
              layout_column_wrap(
                width = 1/2,
                virtualSelectInput("plot_format_Multi", "Format", choices = c("pdf", "svg", "png")),
                numericInput("plot_res_Multi", "Scale", value = 3)
              ),
              downloadButton("download_plot_MultiAbs", "Download Plot", class = "btn-primary w-100")
            )
          )
        ),
        
        # --- 5. RELATIVE ddCq MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Relative ddCq'",
          # Data Selection
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Selection"),
            virtualSelectInput("MockM", "Calibrator (Control)", choices = NULL, placeholder = "Select baseline..."),
            span(class="help-block", "The baseline sample for ddCq."),
            hr(),
            checkboxInput("groupSamplesMulti_dd", "Group Biological Replicates", FALSE),
            hr(),
            virtualSelectInput("SamplePickerDDCtMulti", "Filter Samples", choices = "", multiple = TRUE),
            virtualSelectInput("GenePickerDDCtMulti", "Filter Genes", choices = "", multiple = TRUE)
          ),
          # Data Visualization
          tags$div(
            class = "sidebar-section",
            tags$h6("Data Visualization"),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("PlotDataDDCtMulti", "Data Type", choices = c("Fold Change", "ddCt")),
              virtualSelectInput("PlotTypeDDCtMulti", "Plot Type", choices = c("Bar Chart", "Dot Plot"))
            ),
            layout_column_wrap(
              width = 1/2,
              virtualSelectInput("scaleDDCtMulti", "Y-Axis", choices = c("normal", "log10", "log2")),
              virtualSelectInput("error_typeDDCtMulti", "Error Value", choices = c("SD" = "sd", "SEM" = "sem"))
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("ShowErrorBarDDCtMulti", "Error Bars", FALSE),
              checkboxInput("show_ref_lineDDCtMulti", "Ref. Line", FALSE)
            ),
            layout_column_wrap(
              width = 1/2,
              checkboxInput("show_significanceDDCtMulti", "Significance (*)", FALSE),
              checkboxInput("facet_by_geneDDCtMulti", "Facet by Gene", FALSE)
            )
          ),
          # Appearance
          switchInput("show_appearance_settings_ddCqM", label = "Appearance", value = FALSE, size = "small", onStatus = "success"),
          conditionalPanel(
            condition = "input.show_appearance_settings_ddCqM",
            tags$div(
              class = "sidebar-section",
              tags$h6("Plot Appearance"),
              textInput("plotTitleDDCtMulti", "Plot Title"),
              layout_column_wrap(
                width = 1/2,
                textInput("xlab_customDDCtMulti", "X-Axis Label"),
                textInput("ylab_customDDCtMulti", "Y-Axis Label")
              ),
              virtualSelectInput("colorpickerDDCtMulti", "Color Palette", choices = colorlist),
              layout_column_wrap(
                width = 1/2,
                sliderInput("axis_text_sizeDDCtMulti", "Font Size", min = 6, max = 20, value = 12),
                sliderInput("xTextAngleDDCtMulti", "X-Angle", min = 0, max = 90, value = 0)
              ),
              layout_column_wrap(
                width = 1/2,
                textInput("legendTitleDDCtMulti", "Legend Title"),
                virtualSelectInput("legendPosition_DDCtMulti", "Position", choices = c("right", "none"))
              ),
              orderInput("x_orderDDCtMulti", "Sample Ordering", items = NULL, width = "100%"),
              span(class="help-block", "Drag and drop to reorder.")
            )
          ),
          # Export
          switchInput("show_export_settings_ddCqM", label = "Export Opt.", value = FALSE, size = "small"),
          conditionalPanel(
            condition = "input.show_export_settings_ddCqM",
            tags$div(
              class = "sidebar-section",
              tags$h6("Export Options"),
              layout_column_wrap(
                width = 1/2,
                numericInput("plot_w_DDCtMulti", "Width (in)", value = 7),
                numericInput("plot_h_DDCtMulti", "Height (in)", value = 5)
              ),
              layout_column_wrap(
                width = 1/2,
                virtualSelectInput("plot_format_DDCtMulti", "Format", choices = c("png", "svg", "pdf")),
                numericInput("plot_res_DDCtMulti", "Scale", value = 3)
              ),
              downloadButton("download_plot_MultiRel", "Download Plot", class = "btn-primary w-100")
            )
          )
        ),
        
        # --- 6. STATISTICAL ANALYSIS MULTI ---
        conditionalPanel(
          condition = "input.tabsMulti == 'Statistics'",
          # Analysis Input
          tags$div(
            class = "sidebar-section",
            tags$h6("Analysis Input"),
            virtualSelectInput("limma_inputMulti", "Linear Model Input", choices = c("dCq")),
            tags$div(
              style = "font-size: 0.8rem; padding: 10px; background: #fff3cd; color: #856404; border-radius: 4px; margin-top:10px;",
              tags$b("R Syntax:"), " Samples must start with a letter. Avoid special characters."
            )
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
      navset_card_pill(
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
                card_header(tagList(icon("table"), "Imported Multi-Plate Data")),
                card_body(DTOutput("multiplePlatesData") %>% withSpinner())
              ),
              navset_card_underline(
                nav_panel(
                  title = "Plate View",
                  icon = icon("eye"),
                  card_header(
                    class = "d-flex justify-content-between align-items-center",
                    "Cq per Plate",
                    virtualSelectInput(
                      "PlateSelect", 
                      label = NULL, 
                      choices = "", 
                      placeholder = "Select Plate",
                      width = "180px"
                    )
                  ),
                  plotlyOutput("MultiplotCtCard", height = "500px")
                ),
                nav_panel(
                  title = "Global Distribution",
                  icon = icon("chart-area"),
                  plotlyOutput("MultiplotCtDistrib", height = "500px")
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
                card_header(tagList(icon("filter"), "Well Exclusion Review")),
                card_body(
                  p(class="text-muted small mb-2", "Review flagged wells across all plates. Uncheck to manually re-include a well."),
                  DTOutput("dataMultiplePlatesBadRep") %>% withSpinner()
                )
              ),
              navset_card_underline(
                nav_panel("Outlier Analysis", plotlyOutput("MultiplePlatesFilterPlot", height = "500px")),
                nav_panel("Cq Distribution", plotlyOutput("MultiplePlatesBoxplot", height = "500px")),
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
            
            layout_column_wrap(
              width = 1, # Container über die volle Breite
              fill = FALSE,
              
              layout_column_wrap(
                width = 1/3, # Bestimmt das Verhältnis: 2 Teile Haupt-Card, 1 Teil Methodology
                # Wir nutzen hier fixere Breiten für das 2:1 Verhältnis
                style = "grid-template-columns: 2fr 1fr;", 
                
                # LINKER BEREICH: Die kombinierte Analyse-Card

                    navset_card_underline(
                      # Tab 1: Die Roh-IPCs
                      nav_panel(
                        title = "Inter-Plate Calibrators", 
                        DTOutput("extractedIPCsTable")
                      ),
                      # Tab 2: Der visuelle Check
                      nav_panel(
                        title = "Calibration Check Plot", 
                        plotlyOutput("plotIPCcomparison", height = "550px")
                      ),
                      # Tab 3: Die berechneten Faktoren
                      nav_panel(
                        title = "Calibration Factors",
                        DTOutput("tableCalibrationFactors") %>% withSpinner()
                  )
                ),
                
                # RECHTER BEREICH: Methodology
                card(
                  card_header(tagList(icon("info-circle"), "Methodology")),
                  card_body(
                    uiOutput("ipc_info_text")
                  )
                )
              )
            ),
            
            hr(),
            
            # UNTERER BEREICH: Die resultierenden Daten
            card(
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
                          nav_panel("Ranking Table", DTOutput("ref_ranking_table_multi")),
                          nav_panel("Stability Plot", plotlyOutput("ref_stability_plot_multi"))
                        ),
                    card(
                      fill = FALSE,
                      card_header("Algorithm Info"),
                      card_body(
                        # Hier nutzen wir das gleiche UI-Output wie SP, 
                        # da die Erklärung identisch bleibt
                        uiOutput("algo_description_multi")
                      )
                    )
                  )
        ),
        
        
        # --- TAB 4: RELATIVE dCq ---
        nav_panel(
          title = "Relative dCq",
          icon = icon("chart-column"),
            card(fill = FALSE, 
              card_header(tagList(icon("chart-bar"), "Gene Expression (dCq)")),
              card_body(plotlyOutput("ddctMultiAbsGraph", height = "500px") %>% withSpinner())
            ),
            card(fill = FALSE, 
              card_header(tagList(icon("list"), "Normalized Data Table")),
              card_body(DTOutput("ddctAbsoluteMulti"))
            )
        ),
        
        # --- TAB 5: RELATIVE ddCq ---
        nav_panel(
          title = "Relative ddCq",
          icon = icon("chart-line"),
            card(fill = FALSE, 
              card_header(tagList(icon("chart-line"), "Relative Expression (Fold Change)")),
              card_body(plotlyOutput("ddctRelGraphMulti", height = "500px") %>% withSpinner())
            ),
            card(fill = FALSE, 
              card_header(tagList(icon("clipboard-list"), "ddCq Details")),
              card_body(DTOutput("ddctRelativeMulti"))
            )
        ),
        
        # --- TAB 6: STATISTICAL ANALYSIS ---
        nav_panel(
          title = "Statistics",
          icon = icon("calculator"),
          card(
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
