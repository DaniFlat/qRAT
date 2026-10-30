#' @name qRAT - qPCR Relative Expression Analysis Tool
#' @title Main functions for the qRAT Shiny app
#' @author Daniel Flatschacher
####

####
# Packages
####

#' @import shiny
#' @import ggplot2
#' @importFrom shinyjs runExample alert
#' @importFrom DT dataTableOutput renderDataTable JS
#' @importFrom dplyr last first between summarize
#' @importFrom magrittr extract
#' @importFrom methods show
#' @importFrom methods removeClass
#' @importFrom tools file_ext
#' @importFrom HTqPCR normalizeCtData limmaCtData
#' @importFrom limma makeContrasts
#' @importFrom ggiraph girafe renderGirafe girafeOutput opts_sizing opts_tooltip opts_hover opts_selection opts_zoom opts_toolbar
#' @importFrom shinyWidgets virtualSelectInput
#' @importFrom shinyjqui orderInput


library("shiny")
library("bslib")
library("shinyWidgets")
library("HTqPCR")
library("ddCt")
library("ggiraph")
#library("shinyjs")
#library("ggplot2")
library("scales")
library("data.table")
library("DT")
library("dplyr")
library("waiter")
library("tidyr")
library("magrittr")
library("shinycssloaders")
library("xfun")
library("viridisLite")
library("shinyjqui")
library("ggpubr")
library("rlang")

server <- function(input, output, session) {
  
  waiter_hide()
  MP_files_validated <- reactiveVal(FALSE)

  # Send an input update only when the choices really changed (per session).
  # Without this, every filter change re-sent the choices, reset user selections
  # (reference genes, calibrator, sample order) and re-triggered the downstream chain.
  # Reactive in einem Observer lesen: ein Fehler (z. B. aus einem fehlerhaften Upload) ergibt NULL.
  # Unbehandelte Fehler in Observern beenden in Shiny die ganze Sitzung (grauer Bildschirm);
  # req()-Abbrüche werden unverändert weitergereicht.
  safely <- function(expr) {
    tryCatch(expr, error = function(e) if (inherits(e, "shiny.silent.error")) stop(e) else NULL)
  }

  last_sent_choices <- list()
  send_if_changed <- function(key, choices, update) {
    if (identical(last_sent_choices[[key]], choices)) return(invisible(FALSE))
    last_sent_choices[[key]] <<- choices
    update()
    invisible(TRUE)
  }
  
  
  check_filetype <- function(failed = FALSE) {
    modalDialog(
      title = "File Upload",
      
      if (failed){
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Invalid file type", style = "color: red;"))}
      else {
        removeModal()
      },
      footer = tagList(
        modalButton("Close (Esc)")
      )
    )
  }
  
  
  
  
  
  # --- Hilfsfunktion für den Validierungs-Dialog ---
  SP_validateInputFile <- function(Sample_Column, Well_Column, Gene_Column, Cq_Column) {
    # Zählen, wie viele Checks erfolgreich waren
    success_count <- sum(c(TRUE, Sample_Column, Well_Column, Gene_Column, Cq_Column))
    all_clear <- success_count == 5
    
    # CSS für die Status-Elemente
    status_row <- function(label, is_valid, icon_name) {
      div(
        style = "display: flex; justify-content: space-between; align-items: center; padding: 10px; border-bottom: 1px solid #eee;",
        div(tags$span(icon(icon_name), style = "margin-right: 10px; color: #555;"), tags$b(label)),
        if (is_valid) {
          span(class = "badge bg-success", style = "padding: 8px 12px; border-radius: 20px;", icon("check"), " Found")
        } else {
          span(class = "badge bg-danger", style = "padding: 8px 12px; border-radius: 20px;", icon("xmark"), " Missing")
        }
      )
    }
    
    modalDialog(
      title = tagList(icon("file-circle-check"), " File Validation Report"),
      size = "m",
      easyClose = TRUE,
      
      div(
        style = "padding: 10px;",
        # Header-Status
        if (all_clear) {
          div(class = "alert alert-success", style = "border-radius: 10px; border: none; display: flex; align-items: center;",
              icon("circle-check", class = "fa-2x", style = "margin-right: 15px;"),
              div(h5("Perfect!", style = "margin:0; font-weight: bold;"), p("The file structure is valid and ready for analysis.", style = "margin:0;"))
          )
        } else {
          div(class = "alert alert-warning", style = "border-radius: 10px; border: none; display: flex; align-items: center;",
              icon("triangle-exclamation", class = "fa-2x", style = "margin-right: 15px;"),
              div(h5("Attention Needed", style = "margin:0; font-weight: bold;"), p("Some required columns could not be identified automatically.", style = "margin:0;"))
          )
        },
        
        br(),
        
        # Status Liste
        div(
          style = "background: #f8f9fa; border-radius: 10px; border: 1px solid #e9ecef; overflow: hidden;",
          status_row("File Format (.csv, .txt)", TRUE, "file-code"),
          status_row("Sample Column", Sample_Column, "vial"),
          status_row("Well Column", Well_Column, "border-all"),
          status_row("Gene / Target Column", Gene_Column, "dna"),
          status_row("Cq / Ct Column", Cq_Column, "chart-line")
        ),
        
        if (!all_clear) {
          div(style = "margin-top: 15px; padding: 10px;",
              p(class = "text-muted small", icon("lightbulb"), " Tip: Ensure your column headers match the typical nomenclature (e.g., 'Sample', 'Well', 'Gene', 'Cq') or check for extra headers at the beginning of the file.")
          )
        }
      ),
      
      footer = tagList(
        modalButton("Dismiss", icon = icon("xmark"))
      )
    )
  }
  
  
  
  
  # Marker, ob die aktuellen Dateien bereits validiert wurden
  observeEvent(input$plates, {
    req(input$plates)
    MP_files_validated(FALSE)
    
    ff <- input$plates
    sep <- if (input$SepM == "tab") "\t" else if (input$SepM == "comma") "," else if (input$SepM == "semicolon") ";" else "auto"
    all_results <- list()
    
    for (i in 1:nrow(ff)) {
      # gleiche Einlese-Logik wie die Auswertung (inkl. Geräte-Kopfzeilen); Lesefehler beenden die Sitzung nicht
      rd <- read_plate_file(ff$datapath[i], sep)
      all_results[[i]] <- list(
        plateName = ff$name[i],
        is_sample = is.null(rd$error) && !"Sample" %in% rd$missing,
        is_well   = is.null(rd$error) && !"Well" %in% rd$missing,
        is_gene   = is.null(rd$error) && !"Gene" %in% rd$missing,
        is_cq     = is.null(rd$error) && !"Ct" %in% rd$missing,
        error     = !is.null(rd$error),
        errorMessage = rd$error
      )
    }
    
    MP_validateAllFiles(all_results)
    MP_files_validated(TRUE)
  })
  
  
  
  # --- Hilfsfunktion für den Validierungs-Dialog (Multi Plate) ---
  MP_validateAllFiles <- function(results) {
    
    # Check ob irgendwo ein Fehler vorliegt für das Gesamt-Banner
    any_error <- any(sapply(results, function(x) x$error || !all(x$is_sample, x$is_well, x$is_gene, x$is_cq)))
    
    create_file_box <- function(res) {
      # Wenn die Datei gar nicht lesbar war
      if (isTRUE(res$error)) {
        return(div(
          style = "margin-bottom: 15px; border: 2px solid #d9534f; border-radius: 10px; overflow: hidden; background: #fff;",
          div(style = "background: #d9534f; color: white; padding: 10px 15px; display: flex; justify-content: space-between;", 
              tags$b(res$plateName), icon("circle-xmark")),
          div(style = "padding: 15px; color: #d9534f; background: #fff5f5;", 
              tags$b("Critical Error:"), " This file could not be read at all. Please check if the file is empty, corrupted, or uses a completely different separator.")
        ))
      }
      
      all_clear_file <- all(res$is_sample, res$is_well, res$is_gene, res$is_cq)
      
      # Status Badge Helper für die kompakte Ansicht
      badge_item <- function(label, is_valid) {
        span(
          style = paste0("padding: 4px 10px; border-radius: 15px; font-size: 0.8em; margin-right: 5px; margin-bottom: 5px; display: inline-block; ",
                         if(is_valid) "background: #e8f5e9; color: #2e7d32; border: 1px solid #c8e6c9;" 
                         else "background: #ffegee; color: #c62828; border: 1px solid #ffcdd2;"),
          icon(if(is_valid) "check" else "xmark"), " ", label
        )
      }
      
      div(
        style = paste0("margin-bottom: 15px; border: 1px solid ", if(all_clear_file) "#e0e0e0" else "#f0ad4e", "; border-radius: 10px; overflow: hidden; background: #fff; shadow: 0 2px 4px rgba(0,0,0,0.05);"),
        # Dateiname Header
        div(
          style = paste0("padding: 8px 15px; display: flex; justify-content: space-between; align-items: center; ",
                         if(all_clear_file) "background: #f8f9fa; border-bottom: 1px solid #e0e0e0;" 
                         else "background: #fff8e1; border-bottom: 1px solid #ffe082;"),
          tags$b(res$plateName, style = if(!all_clear_file) "color: #856404;" else "color: #333;"),
          if(all_clear_file) icon("circle-check", style = "color: #28a745;") 
          else icon("triangle-exclamation", style = "color: #f0ad4e;")
        ),
        # Badges Bereich
        div(
          style = "padding: 12px 15px; background: white;",
          badge_item("Sample", res$is_sample),
          badge_item("Well", res$is_well),
          badge_item("Gene", res$is_gene),
          badge_item("Cq/Ct", res$is_cq)
        )
      )
    }
    
    showModal(modalDialog(
      title = tagList(icon("layer-group"), " Multi-Plate Validation Summary"),
      size = "m",
      easyClose = FALSE,
      
      div(
        style = "padding: 5px;",
        # --- Gesamt-Status Banner ---
        if (!any_error) {
          div(class = "alert alert-success", style = "border-radius: 10px; border: none; display: flex; align-items: center; margin-bottom: 20px;",
              icon("circle-check", class = "fa-2x", style = "margin-right: 15px;"),
              div(tags$b("All plates valid!", style = "font-size: 1.1em;"), br(), "All required columns were found in all uploaded files."))
        } else {
          div(class = "alert alert-warning", style = "border-radius: 10px; border: none; display: flex; align-items: center; margin-bottom: 20px;",
              icon("triangle-exclamation", class = "fa-2x", style = "margin-right: 15px;"),
              div(tags$b("Validation issues detected"), br(), "Some files are missing columns or couldn't be read. Please check the details below."))
        },
        
        # --- Liste der Dateien ---
        div(style = "max-height: 400px; overflow-y: auto; padding-right: 5px;",
            lapply(results, create_file_box)),
        
        # --- Tipp Bereich (nur wenn Fehler vorhanden) ---
        if (any_error) {
          div(style = "margin-top: 10px; padding: 10px; background: #fdfdfe; border-radius: 8px; border: 1px dashed #dee2e6;",
              p(class = "text-muted small", style = "margin:0;",
                icon("lightbulb"), " Tip: If all columns are 'Missing', you likely have selected the wrong separator (Comma vs. Semicolon) or the files have extra header lines.")
          )
        }
      ),
      
      footer = tagList(
        modalButton("Confirm & Close", icon = icon("check"))
      )
    ))
  }
  
  
  
  
  
  
  
  ####
  # Data Processing
  ####
  
  #vectors containing typical column names from different qPCR machines; for user-friendly upload of the input data file
  #extend these lists in future versions of qRAT
  sampleColumnList <- c("sample", "sample name", "name")
  wellColumnList <- c("well")
  geneColumnList <- c("gene", "target", "target name", "detector", "primer/probe")
  cqColumnList <- c("ct", "cq")

  # Liest eine Plattendatei ein (Geräte-Kopfzeilen vor der Zeile "Well..." werden übersprungen) und ordnet
  # die Spalten Well/Sample/Gene/Ct zu. Gibt es mehrere passende Spalten (z. B. "Sample" und "Sample Name"),
  # gilt die Reihenfolge der Listen oben. Fehler werden zurückgegeben statt geworfen, damit ein fehlerhafter
  # Upload nicht die Sitzung beendet.
  # Rückgabe: list(data = data.table(Well, Sample, Gene, Ct) oder NULL, missing = fehlende Spalten, error = Lesefehler)
  read_plate_file <- function(path, sep) {
    dt <- tryCatch({
      file.header <- readLines(con = path, n = 100, warn = FALSE)
      n.header <- grep("^Well", file.header) - 1
      if (length(n.header) == 0) n.header <- 0
      fread(path, header = TRUE, sep = sep, skip = n.header[1])
    }, error = function(e) e)
    if (inherits(dt, "error")) return(list(data = NULL, missing = character(0), error = conditionMessage(dt)))

    names_lc <- tolower(names(dt))
    roles <- list(Well = wellColumnList, Sample = sampleColumnList, Gene = geneColumnList, Ct = cqColumnList)
    cols <- vapply(roles, function(cands) {
      hit <- intersect(cands, names_lc)
      if (length(hit)) hit[1] else NA_character_
    }, character(1))
    missing <- names(cols)[is.na(cols)]
    if (length(missing) > 0) return(list(data = NULL, missing = missing, error = NULL))

    out <- dt[, match(cols, names_lc), with = FALSE]
    setnames(out, names(cols))
    # Proben-/Gennamen wie "1" oder "24" sollen Namen bleiben, keine Zahlen
    out[, `:=`(Well = as.character(Well), Sample = as.character(Sample), Gene = as.character(Gene))]
    list(data = out, missing = character(0), error = NULL)
  }

  
  # Single Plate read file
  
  readSinglePlateData <- reactive({

    sep <- switch(input$Sep %||% "auto", tab = "\t", comma = ",", semicolon = ";", "auto")

    req(input$dtfile)
    f <- input$dtfile

    # Check file type - bei falschem Typ abbrechen (vorher wurde trotzdem weitergelesen)
    if (!tolower(file_ext(f$name)) %in% c("csv", "txt", "tsv")) {
      showModal(check_filetype(failed = TRUE))
      return(NULL)
    }

    rd <- read_plate_file(f$datapath, sep)

    ##Validate Input Data File (ein unlesbarer Inhalt zählt wie "alle Spalten fehlen")
    found <- function(role) is.null(rd$error) && !role %in% rd$missing
    showModal(SP_validateInputFile(found("Sample"), found("Well"), found("Gene"), found("Ct")))

    # --- SICHERHEITS-CHECK: ohne alle Pflichtspalten keine Auswertung ---
    if (is.null(rd$data)) return(NULL)
    SP_data <- rd$data
    
    ## generate replicate number for dataview
    SP_data$TempRepNum <- paste(SP_data$Sample, SP_data$Gene)
    SP_data$rp.num <- ave(SP_data$Sample, SP_data$TempRepNum, FUN = seq_along)
    
    # Sicherstellen, dass nur vorhandene Spalten gewählt werden
    SP_data <- subset(SP_data, select = c(Well, Sample, Gene, Ct, rp.num)) %>%
      filter(Well != "")
    
    # ## replace "," with "." for dataview
    # CorrectCt <- SP_data$Ct
    # CorrectCt <- gsub(",", ".", CorrectCt, ignore.case = TRUE)
    # SP_data$Ct <- as.numeric(CorrectCt)
    
    SP_data$Ct <- as.numeric(gsub(",", ".", as.character(SP_data$Ct)))
    
    # Update Inputfield for NTC Selection with all Samplenames
    NTCs <- unique(as.character(SP_data$Sample))
    updateVirtualSelect("NTC_Input", choices = NTCs)
    
    # raw data for plots (spatial plate view; cq distribution)
    plateView <- tryCatch({
      SP_data %>%
        separate(Well,
                 into = c("text", "num"),
                 sep = "(?<=[A-Za-z])(?=[0-9])"
        )
    }, error = function(e) { NULL })
    
    # Falls der Split fehlschlägt, stoppen wir hier, damit app nicht crasht
    if (is.null(plateView)) return(NULL)
    
    return(list(data = SP_data, PV = plateView))
  })
  
  
  
  setData <- reactive({
    req(input$dtfile)
    info <- readSinglePlateData()
    xx <- info$data
    if (is.null(xx)) return(NULL)
    
    experimental_controls <- input$NTC_Input
    
    # --- Step 1: Qaulity ---
    
    # Prepare data
    processed <- xx %>%
      group_by(Sample, Gene) %>%
      mutate(Cq_backup = Ct) %>%
      filter(Sample != "" & Gene != "") %>%
      filter(!Sample %in% experimental_controls)
    
    # --- FILTER Methods ---
    if (input$filterMode == "auto") {
      # AUTOMATIC MODE
      # erst Cq-Grenzen, dann Replikat-Streuung (ein Wert außerhalb der Grenzen soll die SD nicht verfälschen)
      processed <- processed %>%
        mutate(Ct = remove_maxCt(Ct)) %>%
        mutate(Ct = remove_replicates(Ct))
      
    } else {
      # MANUAL MODE
      
      excl_idx <- excluded_rows()
      
      if (length(excl_idx) > 0) {
        # Wir brauchen eine eindeutige ID pro Original-Zeile
        xx_temp <- xx
        xx_temp$OriginalRowID <- 1:nrow(xx_temp)
        
        # Gleicher Filter-Schritt wie oben für 'processed'
        processed_temp <- xx_temp %>%
          filter(Sample != "" & Gene != "") %>%
          filter(!Sample %in% experimental_controls)
        
        # Identifiziere welche dieser Zeilen ausgeschlossen werden sollen
        rows_to_na <- processed_temp$OriginalRowID %in% excl_idx
        
        # Setze Ct auf NA für diese Zeilen
        processed$Ct[rows_to_na] <- NA
      }
    }
    
    processed <- processed %>% ungroup()
    
    # Identifikation der Bad Replicates (Original-Namen behalten!)
    badRep <- processed %>%
      filter(is.na(Ct)) %>%
      rename(Cq = Cq_backup)
    
    # --- Step 2: GROUPING  ---
    grouping_active <- isTRUE(input$groupSamples) && !is.null(confirmedGrouping())
    if (grouping_active) {
      groups <- confirmedGrouping()
      processed$BioSample <- as.character(processed$Sample)   # ursprüngliche (biologische) Probe merken
      new_sample_names <- as.character(processed$Sample)
      
      for (g in groups) {
        grp_name <- g$name
        grp_samples <- g$samples
        
        if (!is.null(grp_name) && nzchar(grp_name) && !is.null(grp_samples)) {
          new_sample_names[new_sample_names %in% grp_samples] <- grp_name
        }
      }
      processed$Sample <- new_sample_names
      processed <- processed %>% mutate(Sample = trimws(Sample)) %>% filter(nzchar(Sample))
    }
    
    # Downstream Updates
    Genes <- unique(processed$Gene)
    send_if_changed("Refs", Genes, function() {
      refs <- grep("^(ACT|UBQ|REF|RF|BACTIN|ACTIN|GAPDH|TUB|TUBULIN).*", Genes, ignore.case = TRUE)
      selected_refs <- if (length(refs) == 0) Genes[length(Genes)] else Genes[refs]
      updateVirtualSelect("Refs", choices = Genes, selected = selected_refs)
    })

    # Auswertungsdaten: bei aktiver Gruppierung geht jede biologische Probe mit dem Mittel ihrer technischen
    # Replikate ein (qPCR.R, aggregate_biological); sonst die technischen Replikate selbst
    analysis <- if (grouping_active) aggregate_biological(processed) else processed

    # qPCR Objekt Erstellung (direkt aus dem Speicher, keine Zwischendatei)
    htset <- NULL; ddset <- NULL; Samples <- NULL
    xdata <- read.qPCRtable(analysis)
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))
      send_if_changed("Mock_Comps", Samples, function() {
        updateVirtualSelect("Mock", choices = Samples, selected = Samples[1])
        updateVirtualSelect("Comps", choices = Samples, selected = Samples)
      })
    }
    
    return(list(
      data = processed,        # technische Replikate (Qualitätskontrolle, Plots der Filterung)
      analysis = analysis,     # Grundlage für dCq/ddCq und Statistik
      data.ht = htset,
      data.ddct = ddset, 
      badReplicates = badRep, 
      CompSamples = Samples
    ))
  })
  
  
  
  # remove samples based on max Ct value, single plates
  remove_maxCt <- function(x, na.rm = TRUE, ...) {
    maxCt <- input$maxCt
    
    y <- x
    y[x < maxCt[1]] <- NA
    y[x > maxCt[2]] <- NA
    y
  }
  
  
  # remove replicates based on threshold, single plates
  # (SD > Grenzwert: der am weitesten vom Mittelwert entfernte Wert wird entfernt, mind. 2 bleiben; qPCR.R)
  remove_replicates <- function(x, na.rm = TRUE, ...) {
    filter_replicates_sd(x, input$repDev)
  }
  
  # remove samples based on max Ct value, multiple plates
  remove_maxCtMulti <- function(x, na.rm = TRUE, ...) {
    maxCt <- input$maxCtMulti
    
    y <- x
    y[x < maxCt[1]] <- NA
    y[x > maxCt[2]] <- NA
    y
  }
  
  
  # remove replicates based on threshold, multiple plates
  remove_replicatesMP <- function(x, na.rm = TRUE, ...) {
    filter_replicates_sd(x, input$repDevMulti)
  }
  
  
  
  # --- Einlese-Funktion (Multi Plate) ---
  readMPData <- reactive({
    req(input$plates)
    
    if (!MP_files_validated()) return(NULL)
    
    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"
    if (input$SepM == "auto") sep <- "auto"
    
    ff <- input$plates
    numfiles <- nrow(ff)
    filelist <- list()
    
    # 1. Dateityp-Check
    for (i in 1:numfiles) {
      ext <- tolower(tools::file_ext(ff$name[i]))
      if (!ext %in% c("csv", "txt", "tsv")) {
        showModal(modalDialog(
          title = "Invalid File Type",
          paste("File", ff$name[i], "is not a supported format (.csv, .txt, .tsv)"),
          easyClose = TRUE
        ))
        return(NULL)
      }
    }
    
    # 2. Dateien einzeln einlesen; Platten ohne Pflichtspalten (oder unlesbare) werden ausgeschlossen,
    #    statt als Platte mit lauter leeren Cq-Werten mitgerechnet zu werden
    skipped <- character(0)
    for (i in 1:numfiles) {
      rd <- read_plate_file(ff$datapath[i], sep)
      if (is.null(rd$data)) {
        reason <- if (!is.null(rd$error)) "not readable" else paste("missing column(s):", paste(sub("^Ct$", "Cq", rd$missing), collapse = ", "))
        skipped <- c(skipped, paste0(ff$name[i], " (", reason, ")"))
      } else {
        filelist[[as.character(i)]] <- rd$data   # Plattennummer = Position der Datei im Upload
      }
    }
    send_if_changed("MP_skipped", skipped, function() {
      if (length(skipped) > 0) {
        showNotification(tagList(tags$b("Plate(s) excluded from the analysis:"), tags$br(), paste(skipped, collapse = "; ")),
                         type = "error", duration = 20)
      }
    })
    if (length(filelist) == 0) return(NULL)

    # 3. Datenverarbeitung
    multiplePlates <- rbindlist(filelist, use.names = TRUE, fill = TRUE, idcol = "PlateNumber")
    multiplePlates$PlateNumber <- as.integer(multiplePlates$PlateNumber)
    
    # Replikate-Nummer generieren
    multiplePlates$TempRepNum <- paste(multiplePlates$Sample, multiplePlates$Gene)
    multiplePlates$rp.num <- ave(multiplePlates$Sample, multiplePlates$TempRepNum, FUN = seq_along)
    
    multiplePlates <- subset(multiplePlates, select = c(PlateNumber, Well, Sample, Gene, Ct, rp.num)) %>%
      filter(Well != "")
    
    # Cq-Werte bereinigen
    multiplePlates$Ct <- as.numeric(gsub(",", ".", as.character(multiplePlates$Ct)))
    
    # Update der UI-Elemente
    ChosenPlate <- unique(as.character(multiplePlates$PlateNumber))
    updateVirtualSelect("PlateSelect", choices = ChosenPlate, selected = ChosenPlate[1])
    
    # Daten für Plate-View vorbereiten
    MultiplePlatesView <- tryCatch({
      multiplePlates %>%
        separate(Well,
                 into = c("text", "num"),
                 sep = "(?<=[A-Za-z])(?=[0-9])"
        )
    }, error = function(e) { NULL }) # Falls Split scheitert, View auf NULL setzen
    
    # Falls View NULL ist (wegen falschem Separator), bricht die Funktion hier sicher ab
    if (is.null(MultiplePlatesView)) return(NULL)
    
    # Update NTC Auswahl
    MP_NTCs <- sort(unique(as.character(multiplePlates$Sample)))
    updateVirtualSelect("NTC_Input_MP", choices = MP_NTCs)
    
    return(list(
      data = multiplePlates, 
      MPV = MultiplePlatesView, 
      all_raw_list = filelist
    ))
  })
  
  
  
  # --- 4. Haupt-Datenverarbeitung (Multi Plate) ---
  multiData <- reactive({
    req(input$plates)
    info <- readMPData()
    multiplePlates <- info$data
    if (is.null(multiplePlates)) return(NULL)
    
    # --- SCHRITT 1: FILTERING & QUALITY ---
    
    # Vorbereitung: ID für manuelles Filtern
    # Wir nutzen eine temporäre Kopie, um die Zeilennummern (OriginalRowID) zu tracken
    temp_df <- multiplePlates
    temp_df$OriginalRowID <- 1:nrow(temp_df)
    
    correctRep <- temp_df %>%
      group_by(PlateNumber, Sample, Gene) %>%
      mutate(Cq_backup = Ct) %>%
      filter(Sample != "" & Gene != "") %>%
      filter(!Sample %in% input$NTC_Input_MP)
    
    # --- FILTER Methods---
    if (input$filterModeMP == "auto") {
      # AUTOMATIC MODE
      correctRep <- correctRep %>%
        mutate(Ct = remove_maxCtMulti(Ct)) %>%
        mutate(Ct = remove_replicatesMP(Ct))
      
    } else {
      # MANUAL MODE
      excl_idx <- excluded_rows_MP()
      
      if (length(excl_idx) > 0) {
        # Identifiziere Zeilen basierend auf der Original-ID
        rows_to_na <- correctRep$OriginalRowID %in% excl_idx
        
        # Setze Ct auf NA
        correctRep$Ct[rows_to_na] <- NA
      }
    }
    
    correctRep <- correctRep %>% ungroup()
    
    # Bad Replicates definieren (wo Ct NA ist)
    badRep <- correctRep %>% 
      filter(is.na(Ct)) %>% 
      rename(Cq = Cq_backup)
    
    # --- SCHRITT 2: SAMPLE GROUPING ---
    if (isTRUE(input$groupSamplesMulti) && !is.null(confirmedGroupingMulti())) {
      groups <- confirmedGroupingMulti()
      correctRep$BioSample <- as.character(correctRep$Sample)   # ursprüngliche (biologische) Probe merken
      new_names <- as.character(correctRep$Sample)
      
      for (g in groups) {
        grp_name <- g$name
        grp_samples <- g$samples
        if (!is.null(grp_name) && !is.null(grp_samples)) {
          new_names[new_names %in% grp_samples] <- grp_name
        }
      }
      correctRep$Sample <- new_names
    }
    
    # Update IPC Selection
    SamplesWithIPCs <- sort(unique(as.character(correctRep$Sample)))
    # isolate(): a change of the IPC selection must not re-run the whole filtering step
    send_if_changed("IPC", SamplesWithIPCs, function() {
      updateVirtualSelect("IPC", choices = SamplesWithIPCs, selected = isolate(input$IPC))
    })
    
    # Spalte OriginalRowID entfernen, falls sie später stört
    correctRep$OriginalRowID <- NULL
    
    return(list(data = correctRep, badReplicates = badRep))
  })
  
  # --- multiData 1b: Trennung der IPCs von den Haupt-Proben (unkalibrierte Cq-Werte) ---
  multiDataNoIPCs <- reactive({
    req(input$plates)

    info <- multiData()
    if (is.null(info)) return(NULL)

    info2 <- info$data
    if (is.null(info2)) return(NULL)

    # Falls IPC Kalibrierung aktiv ist, IPC Probe aus dem dCt/ddCt Datensatz entfernen
    if (isTRUE(input$Id027)) {
      req(input$IPC)
      ipcs <- input$IPC
      info2 <- info2[!(info2$Sample %in% ipcs), ]
    }

    info2
  })

  # --- multiData 2: qPCR-Objekte für dCq/ddCq und limma ---
  multiData2 <- reactive({
    info2 <- multiDataNoIPCs()
    if (is.null(info2)) return(NULL)

    # Bei aktiver IPC-Kalibrierung werden die kalibrierten Cq-Werte ausgewertet
    calibrated <- if (isTRUE(input$Id027)) interPlateCalibration() else info2
    # bei aktiver Gruppierung zählt jede biologische Probe einmal (Mittel ihrer technischen Replikate, über alle Platten)
    analysis_data <- if ("BioSample" %in% names(calibrated)) aggregate_biological(calibrated) else calibrated

    # Gene einlesen und Housekeeper vorschlagen
    Genes <- unique(info2$Gene)
    send_if_changed("RefsM", Genes, function() {
      refs <- grep("^(ACT|UBQ|REF|RF|BACTIN|ACTIN|GAPDH|TUB|TUBULIN).*", Genes, ignore.case = TRUE)
      selected_refs <- if (length(refs) == 0) Genes[length(Genes)] else Genes[refs]
      updateVirtualSelect("RefsM", choices = Genes, selected = selected_refs)
    })

    # Konvertierung in das interne qPCR-Format (direkt aus dem Speicher, keine Zwischendatei)
    htset <- NULL; ddset <- NULL; Samples <- NULL
    xdata <- read.qPCRtableMulti(analysis_data)

    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))

      send_if_changed("MockM_CompsM", Samples, function() {
        updateVirtualSelect("MockM", choices = Samples, selected = Samples[1])
        updateVirtualSelect("CompsM", choices = Samples, selected = Samples)
      })
    }
    
    return(list(
      SamplesNoIPCs = info2,
      SamplesCalibrated = calibrated,   # je Well (mit Plattennummer), z. B. für den Reference Finder
      SamplesAnalysed = analysis_data,
      xdata = xdata,
      data.ht = htset, 
      data.ddct = ddset, 
      CompSamples = Samples
    ))
  })
  
  
  
  
  # UI-Auswahl der p-Wert-Korrektur ("BH", "Holm", "Bonferroni") auf die limma-Methodennamen abbilden.
  # Ohne Auswahl gilt BH (erster Eintrag im Dropdown).
  limma_adjust_method <- function(choice) {
    if (identical(choice, "Holm")) {
      "holm"
    } else if (identical(choice, "Bonferroni")) {
      "bonferroni"
    } else {
      "BH"
    }
  }

  ## HTqPCR Single Plate
  resultHTqPCR <- reactive({
    info <- setData()

    if (is.null(info)) {
      return(NULL)
    }

    # "dCq" ist die einzige Option für input$limma_input
    dataset <- info$data.ht

    refs <- input$Refs
    multiComps <- listOfCompsMultiple()

    if (input$compType == "Single Comparison") {
      comp.type <- 1
      comps <- input$Comps
    } else {
      comp.type <- 2
      comps <- multiComps
    }

    adjustMethod <- limma_adjust_method(input$adjustMethod)

    res <- calHTqPCR(dataset, refs, comps, comp.type, adjustMethod)
    return(res)
  })
  
  ## HTqPCR Multiple Plates
  resultHTqPCRMulti <- reactive({
    info <- multiData2()
    if (is.null(info)) {
      return(NULL)
    }
    
    # "dCq" ist die einzige Option für input$limma_inputMulti
    dataset <- info$data.ht

    refs <- input$RefsM
    multiComps <- listOfCompsMultipleM()
    
    if (input$compTypeM == "Single Comparison") {
      comp.type <- 1
      comps <- input$CompsM
    } else {
      comp.type <- 2
      comps <- multiComps
    }
    
    adjustMethod <- limma_adjust_method(input$adjustMethodMulti)
    
    res <- calHTqPCR(dataset, refs, comps, comp.type, adjustMethod)
    return(res)
  })
  
  ## Limma (HTqPCR) Single Plate
  calLimma <- reactive({
    info <- resultHTqPCR()
    if (is.null(info)) {
      return(NULL)
    }
    if (is.null(info$result)) {
      return(NULL)
    }
    cts <- info$contrast
    cns <- names(cts)
    df <- info$df
    ## data.frame(cns)
    res <- NULL
    for (i in 1:length(cts)) {
      aa <- info$result[[i]]
      aa$Label <- rep(cns[i], nrow(aa))
      aa <- aa[, c("Label", "genes", "FC", "t.test", "p.value", "adj.p.value")]
      res <- rbind(res, aa)
    }
    res <- res[order(res$Label, res$genes), ]
    res[!res$genes %in% input$Refs, ]
  })
  
  ## Limma (HTqPCR) Multiple Plates
  calLimmaMulti <- reactive({
    info <- resultHTqPCRMulti()
    if (is.null(info)) {
      return(NULL)
    }
    if (is.null(info$result)) {
      return(NULL)
    }
    cts <- info$contrast
    cns <- names(cts)
    ## data.frame(cns)
    res <- NULL
    for (i in 1:length(cts)) {
      aa <- info$result[[i]]
      aa$Label <- rep(cns[i], nrow(aa))
      aa <- aa[, c("Label", "genes", "FC", "t.test", "p.value", "adj.p.value")]
      res <- rbind(res, aa)
    }
    res <- res[order(res$Label, res$genes), ]
    res[!res$genes %in% input$RefsM, ]
  })
  
  
  
  
  ## ddCt
  # Anzahl gültiger Cq-Replikate des Zielgens je Sample/Gene (gleiche Regel wie read.qPCRtable: Cq > 40 gilt als NA)
  add_replicate_n <- function(result, data, na.value = 40) {
    counts <- data %>% filter(!is.na(Ct), Ct <= na.value) %>% count(Sample, Gene)
    key <- function(s, g) paste(as.character(s), as.character(g), sep = "\r")
    n <- counts$n[match(key(result$Sample, result$Gene), key(counts$Sample, counts$Gene))]
    result$n <- ifelse(is.na(n), 0L, n)
    result
  }

  # Hinweis auf Werte ohne SD (Ziel- oder Referenzgen mit < 2 gültigen Replikaten). Früher wurden diese
  # Zeilen stillschweigend entfernt; jetzt bleiben sie ohne Fehlerbalken in Tabelle und Plot.
  notify_missing_sd <- function(key, expr.abs, expr.rel) {
    missing <- unique(c(
      paste(expr.abs$Sample, expr.abs$Gene, sep = " / ")[is.na(expr.abs$dCt.sd)],
      paste(expr.rel$Sample, expr.rel$Gene, sep = " / ")[is.na(expr.rel$ddCt.sd)]
    ))
    send_if_changed(key, missing, function() {
      if (length(missing) > 0) {
        shown <- if (length(missing) > 10) c(head(missing, 10), "...") else missing
        showNotification(
          tagList(
            tags$b(paste0("No standard deviation for ", length(missing), " sample/gene combination(s):")),
            tags$br(), paste(shown, collapse = ", "), tags$br(),
            "The target or reference gene has fewer than 2 valid replicates. These values are shown without error bars."
          ),
          type = "warning", duration = 20
        )
      }
    })
  }

  resultDdct <- reactive({
    info <- setData()
    if (is.null(info)) {
      return(NULL)
    }
    ddct.raw <- info$data.ddct
    if (is.null(ddct.raw)) {
      return(NULL)
    }
    mock <- input$Mock
    refs <- input$Refs
    
    result.dd <- ddCtExpression(ddct.raw, calibrationSample = mock, housekeepingGene = refs)
    result.dd <- elist(result.dd)
    
    expr.rel <- result.dd[, c("Sample", "Detector", "exprs", "level.err", "ddCt", "ddCt.error")]
    colnames(expr.rel) <- c("Sample", "Gene", "FC", "FC.sd", "ddCt", "ddCt.sd")
    expr.rel <- expr.rel[!expr.rel$Gene %in% refs, ] %>%
      filter(Sample != "") %>%
      filter(Gene != "") %>%
      filter(!is.na(FC)) # Zeilen ohne SD (n < 2) bleiben erhalten
    
    expr.abs <- result.dd[, c("Sample", "Detector", "dCt", "dCt.error")] %>%
      mutate(dCt2 = dCt, dCt3 = dCt, dCt.error2 = dCt.error) # duplicate dCt results column, for transforming one of the duplicates into expressions
    colnames(expr.abs) <- c("Sample", "Gene", "dCt", "dCt.sd", "-dCq", "expr", "expr.sd")
    aa <- expr.abs$expr
    bb <- expr.abs$expr.sd
    expr.abs$expr <- 2^(-aa) #* 1000
    expr.abs$expr.sd <- abs(expr.abs$expr - 2^(-(aa + bb))) #* 1000)
    expr.abs <- expr.abs[!expr.abs$Gene %in% refs, ] %>%
      filter(Sample != "") %>%
      filter(Gene != "") %>%
      filter(!is.na(expr)) # Zeilen ohne SD (n < 2) bleiben erhalten
    expr.abs$"-dCq" <- expr.abs$"-dCq" * -1

    # Spalte n: gültige Replikate des Zielgens (Basis für SEM und Hinweis bei fehlender SD)
    expr.abs <- add_replicate_n(expr.abs, info$analysis)
    expr.rel <- add_replicate_n(expr.rel, info$analysis)
    notify_missing_sd("missing_sd_single", expr.abs, expr.rel)

    SamplesdCt <- unique(as.character(expr.abs$Sample))
    SamplesddCt <- unique(as.character(expr.rel$Sample))
    GenesdCt <- unique(as.character(expr.abs$Gene))
    GenesddCt <- unique(as.character(expr.rel$Gene))
    
    # Update Inputfields (Plot Settings and Input Settings) - only when the lists change,
    # so user selections and the custom sample order survive a change of Refs/Mock
    send_if_changed("SamplePicker", SamplesdCt, function() {
      updateVirtualSelect("SamplePicker", choices = SamplesdCt, selected = SamplesdCt)
      updateOrderInput(session, "x_order", items = SamplesdCt)
    })
    send_if_changed("SamplePickerDDCt", SamplesddCt, function() {
      updateVirtualSelect("SamplePickerDDCt", choices = SamplesddCt, selected = SamplesddCt)
      updateOrderInput(session, "x_orderDDCt", items = SamplesddCt)
    })
    send_if_changed("GenePicker", GenesdCt, function() {
      updateVirtualSelect("GenePicker", choices = GenesdCt, selected = GenesdCt)
    })
    send_if_changed("GenePickerDDCt", GenesddCt, function() {
      updateVirtualSelect("GenePickerDDCt", choices = GenesddCt, selected = GenesddCt)
    })

    # The former ddCt -> "ddCtsingle.csv" -> qPCRset round trip (stat_set) was never used; removed
    list(relative = expr.rel, absolute = expr.abs)
  })
  
  
  
  
  ## ddCt MultiPlate
  resultDdctMulti <- reactive({
    req(input$plates)
    
    info <- multiData2()
    if (is.null(info)) {
      return(NULL)
    }
    ddct.raw <- info$data.ddct
    if (is.null(ddct.raw)) {
      return(NULL)
    }
    mock <- input$MockM
    refs <- input$RefsM
    
    
    result.dd <- ddCtExpression(ddct.raw, calibrationSample = mock, housekeepingGene = refs)
    result.dd <- elist(result.dd)
    resultElist <- result.dd
    expr.rel <- result.dd[, c("Sample", "Detector", "exprs", "level.err", "ddCt", "ddCt.error")]
    colnames(expr.rel) <- c("Sample", "Gene", "FC", "FC.sd", "ddCt", "ddCt.sd")
    
    expr.rel <- expr.rel[!expr.rel$Gene %in% refs, ] %>%
      filter(Sample != "") %>%
      filter(Gene != "") %>%
      filter(!is.na(FC)) # Zeilen ohne SD (n < 2) bleiben erhalten
    
    expr.abs <- result.dd[, c("Sample", "Detector", "dCt", "dCt.error")] %>%
      mutate(dCt2 = dCt, dCt3 = dCt, dCt.error2 = dCt.error) # duplicate dCt results column, for transforming one of the duplicates into expressions
    colnames(expr.abs) <- c("Sample", "Gene", "dCt", "dCt.sd", "-dCq", "expr", "expr.sd")
    aa <- expr.abs$expr
    bb <- expr.abs$expr.sd
    expr.abs$expr <- 2^(-aa) #* 1000
    expr.abs$expr.sd <- abs(expr.abs$expr - 2^(-(aa + bb))) #* 1000)
    
    expr.abs <- expr.abs[!expr.abs$Gene %in% refs, ] %>%
      filter(Sample != "") %>%
      filter(Gene != "") %>%
      filter(!is.na(expr)) # Zeilen ohne SD (n < 2) bleiben erhalten
    
    expr.abs$"-dCq" <- expr.abs$"-dCq" * -1

    # Spalte n: gültige Replikate des Zielgens über alle Platten (bei IPC-Kalibrierung: kalibrierte Daten)
    expr.abs <- add_replicate_n(expr.abs, info$SamplesAnalysed)
    expr.rel <- add_replicate_n(expr.rel, info$SamplesAnalysed)
    notify_missing_sd("missing_sd_multi", expr.abs, expr.rel)

    SamplesdCtMulti <- unique(as.character(expr.abs$Sample))
    SamplesddCtMulti <- unique(as.character(expr.rel$Sample))
    
    GenesdCtMulti <- unique(as.character(expr.abs$Gene))
    GenesddCtMulti <- unique(as.character(expr.rel$Gene))
    
    # Update Inputfields (Plot Settings and Input Settings) - only when the lists change
    send_if_changed("SamplePickerMulti", SamplesdCtMulti, function() {
      updateVirtualSelect("SamplePickerMulti", choices = SamplesdCtMulti, selected = SamplesdCtMulti)
      updateOrderInput(session, "x_orderMulti", items = SamplesdCtMulti)
    })
    send_if_changed("SamplePickerDDCtMulti", SamplesddCtMulti, function() {
      updateVirtualSelect("SamplePickerDDCtMulti", choices = SamplesddCtMulti, selected = SamplesddCtMulti)
      updateOrderInput(session, "x_orderDDCtMulti", items = SamplesddCtMulti)
    })
    send_if_changed("GenePickerMulti", GenesdCtMulti, function() {
      updateVirtualSelect("GenePickerMulti", choices = GenesdCtMulti, selected = GenesdCtMulti)
    })
    send_if_changed("GenePickerDDCtMulti", GenesddCtMulti, function() {
      updateVirtualSelect("GenePickerDDCtMulti", choices = GenesddCtMulti, selected = GenesddCtMulti)
    })
    
    ## updateRadioButtons(session, 'geneNameRel', choices=targets, selected=targets[1])
    ## updateCheckboxGroupInput(session, 'geneNameAbs', choices=targets, selected=targets)
    list(multiRelative = expr.rel, multiAbsolute = expr.abs, elist = resultElist)
  })
  
  
  
  
  extractIPCs <- reactive({
    req(input$Id027, input$IPC)
    if (input$IPC == "") return(NULL)
    
    info <- multiData()
    req(info$data)
    info2 <- info$data
    ipcs <- input$IPC
    
    # 1. IPCs extrahieren
    extractedIPCs <- info2[info2$Sample %in% ipcs, ]
    extractedIPCs <- na.omit(extractedIPCs)
    req(nrow(extractedIPCs) > 0)
    
    # 2. Globalen Mittelwert aller IPCs berechnen (arithmetisch, da Ct log ist)
    # Das arithmetische Mittel von log-Werten entspricht dem geom. Mittel der Rohwerte.
    global_ipc_mean <- mean(extractedIPCs$Ct, na.rm = TRUE)
    
    # 3. Kalibrierungs-Faktoren berechnen
    # Wir berechnen pro Platte: Wie weit weicht der IPC-Schnitt vom globalen Schnitt ab?
    CalibrationFactors <- extractedIPCs %>%
      group_by(PlateNumber) %>%
      summarize(
        Plate_IPC_Mean = mean(Ct, na.rm = TRUE),
        # Der 'Shift' oder 'Mean' Offset:
        Mean = Plate_IPC_Mean - global_ipc_mean 
      )
    
    return(list(
      extrIPCs = extractedIPCs, 
      CalibrationFactors = CalibrationFactors,
      GlobalMean = global_ipc_mean
    ))
  })
  
  
  
  
  
  interPlateCalibration <- reactive({
    # basiert auf den unkalibrierten Daten (nicht auf multiData2, das selbst die kalibrierten Werte nutzt)
    x <- multiDataNoIPCs()
    req(x)

    if (isTRUE(input$Id027)) {
      calF <- extractIPCs()
      req(calF$CalibrationFactors)

      # Platten ohne IPC-Werte erhalten keinen Kalibrierungsfaktor -> Cq wird dort NA
      missing_plates <- setdiff(unique(x$PlateNumber), calF$CalibrationFactors$PlateNumber)
      send_if_changed("IPC_missing_plates", missing_plates, function() {
        if (length(missing_plates) > 0) {
          showNotification(
            paste0("No IPC values found on plate(s) ", paste(missing_plates, collapse = ", "),
                   ": these plates cannot be calibrated, their Cq values are set to NA."),
            type = "warning", duration = 15
          )
        }
      })

      x %>%
        left_join(calF$CalibrationFactors, by = "PlateNumber") %>%
        # Korrektur: Ct_neu = Ct_alt - Offset
        mutate(Ct = Ct - Mean) %>% 
        select(-Mean, -Plate_IPC_Mean)
    } else {
      return(x)
    }
  })
  
  
  # Daten für den Vergleich
  interPlateCalibration_Plot <- reactive({
    calibrated_df <- interPlateCalibration()
    req(calibrated_df)
    
    noncalibrated_df <- multiDataNoIPCs()
    req(noncalibrated_df)
    
    # Hilfsspalte für die Auswahl erstellen
    cal <- calibrated_df %>% unite("Sample_Gene", Sample:Gene, remove = FALSE)
    non <- noncalibrated_df %>% unite("Sample_Gene", Sample:Gene, remove = FALSE)
    
    return(list(calibrated = cal, noncalibrated = non))
  })
  
  # Separater Observer für das Update des Dropdowns
  observe({
    info <- safely(interPlateCalibration_Plot())
    req(info)
    
    choices_list <- unique(info$calibrated$Sample_Gene)
    updateVirtualSelect("SamplePickerIPCcomparison", choices = choices_list)
  })
  
  
  
  # table options
  table_options <- function() {
    list(
      dom = "Brtip",
      # Bfrtip
      pageLength = 10,
      row_number = FALSE,
      buttons =
        list(
          "copy", list(
            extend = "collection",
            buttons = c("csv", "excel", "pdf"),
            text = "Export"
          ),
          list(
            extend = "",
            text = "Show All",
            action = JS(
              "function ( e, dt, node, config ) {
          dt.page.len(-1);
          dt.ajax.reload();}"
            )
          ),
          list(
            extend = "",
            text = "Show Less",
            action = JS(
              "function ( e, dt, node, config ) {
          dt.page.len(10);
          dt.ajax.reload();}"
            )
          )
        ),
      deferRender = TRUE,
      lengthMenu = list(c(10, 20, -1), c("10", "20", "All")),
      digits = 4,
      editable = FALSE,
      scroller = TRUE,
      autoWidth = TRUE,
      processing = FALSE,
      lengthChange = FALSE,
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#325D88', 'color': '#fff'});",
        "}"
      )
    )
  }
  
  
  ## Results Output
  
  ## SINGLE PLATE OUTPUT
  
  ## Raw Data Table Single Plate Original Data
  
  
  output$dataSinglePlate <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- readSinglePlateData()
    req(info)
    
    info$data  %>% rename(Replicate = rp.num)
    
  })
  
  # --- Plot-Bausteine für Rohdaten & Qualität (ggiraph; Single und Multi Plate nutzen dieselben Funktionen) ---

  qc_theme <- function() {
    theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(face = "bold"))
  }

  # Plate View: Heatmap der Cq-Werte je Well. data_id = Zeilennummer in den Rohdaten, damit
  # die per Klick ausgewählten Wells direkt als manuelle Ausschlüsse übernommen werden können.
  plate_view_plot <- function(pv, row_id, selected = character(0)) {
    pv <- as.data.frame(pv)
    pv$row_id <- as.character(row_id)
    pv$col <- as.integer(pv$num)
    pv$row_letter <- factor(pv$text, levels = rev(sort(unique(pv$text))))
    pv$tooltip <- paste0("<b>Well ", pv$text, pv$num, "</b><br>Sample: ", pv$Sample, "<br>Gene: ", pv$Gene,
                         "<br>Cq: ", ifelse(is.na(pv$Ct), "NA", round(pv$Ct, 2)))
    p <- ggplot(pv, aes(x = col, y = row_letter, fill = Ct, tooltip = tooltip, data_id = row_id)) +
      geom_tile_interactive(color = "white", linewidth = 0.6) +
      scale_fill_distiller(palette = "RdYlBu", direction = 1, name = "Cq", na.value = "grey85") +
      scale_x_continuous(position = "top", breaks = sort(unique(pv$col)), expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_fixed() +
      labs(x = NULL, y = NULL, title = "Plate View") +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text = element_text(size = 12, colour = "grey20"))
    ratio <- length(unique(pv$text)) / length(unique(pv$col))
    qrat_girafe(p, width_svg = 6, height_svg = 6 * ratio + 1.2,
                selection = opts_selection(type = "multiple", selected = selected,
                                           css = "stroke:#000000;stroke-width:3px;opacity:0.55;"))
  }

  # Cq-Verteilung je Probe; Hover hebt die Kurve einer Probe hervor
  cq_density_plot <- function(df) {
    df <- as.data.frame(df) %>% filter(!is.na(Ct))
    n <- length(unique(df$Sample))
    # Set1 hat nur 9 Farben; bei mehr Proben wird interpoliert (vorher blieben diese Kurven farblos)
    cols <- if (n <= 9) RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)] else colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
    p <- ggplot(df, aes(x = Ct, colour = Sample)) +
      geom_density_interactive(aes(tooltip = Sample, data_id = Sample), linewidth = 0.8) +
      geom_rug() +
      scale_colour_manual(values = cols) +
      labs(title = "Cq Distribution", x = "Cq", y = "Density", colour = "Samples") +
      qc_theme() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    qrat_girafe(p, width_svg = 6, height_svg = 4.5)
  }

  # Anzahl verwendeter (OK) und verworfener (unreliable) Replikate je Probe
  replicate_status_plot <- function(bad, data) {
    filtered <- bad %>% count(Sample)
    unfiltered <- data %>% filter(!is.na(Ct)) %>% count(Sample)
    together <- bind_rows(filtered, unfiltered, .id = "status") %>%
      mutate(status = factor(ifelse(status == 1, "unreliable", "OK"), levels = c("OK", "unreliable")))
    p <- ggplot(together, aes(x = Sample, y = n, fill = status,
                              tooltip = paste0("<b>", Sample, "</b><br>", status, ": ", n), data_id = paste(Sample, status))) +
      geom_col_interactive(position = position_stack(reverse = TRUE), color = "black", linewidth = 0.5) +
      scale_fill_brewer(palette = "Paired", name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(title = "filtered/unfiltered Replicates [n]", x = NULL, y = "Count") +
      qc_theme()
    qrat_girafe(p, width_svg = 6, height_svg = 4.5)
  }

  # Cq-Werte je Probe und Gen (Box + einzelne Replikate)
  replicate_boxplot <- function(data) {
    df <- as.data.frame(data) %>% filter(!is.na(Ct))
    plate <- if ("PlateNumber" %in% names(df)) paste0("Plate ", df$PlateNumber, ", ") else ""
    df$tooltip <- paste0("<b>", df$Sample, " / ", df$Gene, "</b><br>", plate, "Well ", df$Well, "<br>Cq: ", round(df$Ct, 2))
    df$point_id <- paste0("pt", seq_len(nrow(df)))
    genes <- unique(as.character(df$Gene))
    cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(max(length(genes), 2))[seq_along(genes)]
    p <- ggplot(df, aes(x = Sample, y = Ct, fill = Gene)) +
      geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.6, color = "grey25") +
      geom_point_interactive(aes(tooltip = tooltip, data_id = point_id),
                             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1),
                             shape = 21, size = 1.8, color = "black") +
      scale_fill_manual(values = setNames(cols, genes)) +
      labs(title = "Variation within Samples", x = NULL, y = "Cq", fill = "Genes") +
      qc_theme()
    qrat_girafe(p, width_svg = 6, height_svg = 4.5)
  }


  ## Raw Data Table Plate Plot Spatial Original Data

  output$plotCtCard <- renderGirafe({
    info <- readSinglePlateData()
    req(info)
    # im manuellen Modus sind die ausgeschlossenen Wells vorausgewählt
    plate_view_plot(info$PV, row_id = seq_len(nrow(info$PV)), selected = as.character(excluded_rows()))
  })

  # Ausgewählte Wells als manuelle Ausschlüsse übernehmen
  observeEvent(input$apply_plate_selection, {
    selected <- sort(as.integer(input$plotCtCard_selected))
    excluded_rows(selected)
    updateRadioButtons(session, "filterMode", selected = "manual")
    showNotification(paste(length(selected), "wells excluded (Manual Filtering)."), type = "message")
  })

  ## Raw Data Table Cq distribution Original Data

  output$plotCtDistrib <- renderGirafe({
    info <- readSinglePlateData()
    req(info)
    cq_density_plot(info$PV)
  })
  
  
  ## Single Plate Show bad replicates Table
  
  output$dataSinglePlateBadRep <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- setData()
    req(info)
    data <- info$badReplicates %>%
      select(Well, Sample, Gene, rp.num, Cq) %>%
      rename(Replicate = rp.num)
    
    if (is.null(info)) NULL else data
  })
  
  
  ## Single Plate Filter Plot (how many replicates are neglected per sample)
  
  output$SinglePlateFilterPlot <- renderGirafe({
    info <- setData()
    req(info)
    replicate_status_plot(info$badReplicates, info$data)
  })



  ## Single Plate Boxplot after filtering

  output$SinglePlateBoxplot <- renderGirafe({
    info <- setData()
    req(info)
    replicate_boxplot(info$data)
  })
  
  ## dCt Expression Data Table Single Plate
  
  output$ddctAbsolute <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdct()
    if (is.null(info)) {
      return(NULL)
    }
    
    info$absolute %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      rename(dCq = dCt, dCq.sd = dCt.sd, RQ = expr, RQ.sd = expr.sd)
  })
  

  
  # --- Zentrale qPCR Plot-Funktion ---
  
  render_qPCR_plot <- function(df, stats_df, PlotDataCol, PlotDataError, yTitle, input_type = "dCq", export_mode = FALSE) {
    
    # Hilfsfunktion für verschiedene UI-Suffixe
    get_input <- function(id) {
      # Mapping der Typen auf die UI-Suffixe
      suffix <- switch(input_type,
                       "dCq"      = "",
                       "ddCq"     = "DDCt",
                       "MultiAbs" = "Multi",
                       "MultiRel" = "DDCtMulti",
                       # Suffixe, unter denen der Download-Handler die Multi-Plate-Plots exportiert
                       "Multi"     = "Multi",
                       "DDCtMulti" = "DDCtMulti",
                       "" # Default
      )
      
      # Suche erst mit Suffix, dann ohne
      val <- input[[paste0(id, suffix)]]
      if (is.null(val)) val <- input[[id]]
      return(val)
    }
    
    # 1. Parameter laden
    PlotType       <- get_input("PlotType") %||% "Bar Chart"
    colorPick      <- get_input("colorpicker") %||% "Paired"
    axis_text_size <- get_input("axis_text_size") %||% 10
    xTextAngle     <- get_input("xTextAngle") %||% 0
    plotTitle      <- get_input("plotTitle") %||% ""
    legendTitle    <- get_input("legendTitle") %||% "Gene"
    legendPosition <- get_input("legendPosition") %||% "right"
    scalePick      <- get_input("scale") %||% "normal"
    # (Exportformat und -größe liest der Download-Handler aus plot_format_*/plot_w_*/plot_h_*/plot_res_*)

    # 2. Achsen Titel
    xTitle <- if(nzchar(get_input("xlab_custom") %||% "")) get_input("xlab_custom") else "Sample"
    final_yTitle <- if(nzchar(get_input("ylab_custom") %||% "")) get_input("ylab_custom") else yTitle
    
    # 3. Sample Ordering
    x_order <- get_input("x_order")
    if (!is.null(x_order) && length(x_order) > 0) {
      existing_samples <- intersect(x_order, unique(df$Sample))
      if (length(existing_samples) > 0) {
        df$Sample <- factor(df$Sample, levels = existing_samples)
      }
    } else {
      df$Sample <- as.factor(df$Sample)
    }
    
    # 4. Farb-Management
    numberOfGenes <- n_distinct(df$Gene)
    GeneNames     <- as.character(unique(df$Gene))
    
    if (colorPick == "viridis") {
      colourpalette <- viridisLite::viridis(numberOfGenes)
    } else {
      min_colors <- 3
      max_colors <- RColorBrewer::brewer.pal.info[colorPick, "maxcolors"] %||% 8
      pal_func <- colorRampPalette(RColorBrewer::brewer.pal(min(max(min_colors, numberOfGenes), max_colors), colorPick))
      colourpalette <- pal_func(numberOfGenes)
    }
    names(colourpalette) <- GeneNames
    
    # 5. Tooltip-Vorbereitung
    val_vec <- as.numeric(df[[PlotDataCol]])
    err_vec <- as.numeric(df[[PlotDataError]])
    
    df$tooltip_text <- paste0(
      "<b>Sample:</b> ", df$Sample,
      "<br><b>Gene:</b> ", df$Gene,
      "<br><b>Value:</b> ", round(val_vec, 3),
      "<br><b>Error:</b> ", ifelse(is.na(err_vec), "&ndash; (no SD)", round(err_vec, 3)),
      if (!is.null(df$n)) paste0("<br><b>n:</b> ", df$n) else ""
    )
    df$hover_id <- paste(df$Sample, df$Gene, sep = " / ")

    # 6. Plot-Basis erstellen (interaktive ggiraph-Geoms; beim Export werden sie wie normale Geoms gezeichnet)
    p <- ggplot(df, aes(x = Sample, y = !!sym(PlotDataCol), fill = Gene, color = Gene, group = Gene,
                        tooltip = tooltip_text, data_id = hover_id))

    if (PlotType == "Bar Chart") {
      p <- p + geom_bar_interactive(stat = "identity", position = position_dodge(0.9), width = 0.7)
    } else if (PlotType == "Boxplot") {
      p <- p + geom_boxplot_interactive(position = position_dodge(0.9))
    } else {
      p <- p + geom_point_interactive(position = position_dodge(0.9), size = 3)
    }
    
    # 7. Error Bars
    if (isTRUE(get_input("ShowErrorBar"))) {
      p <- p + geom_errorbar(
        aes(
          ymin = !!sym(PlotDataCol) - !!sym(PlotDataError), 
          ymax = !!sym(PlotDataCol) + !!sym(PlotDataError)
        ),
        width = 0.25,
        color = "black",
        position = position_dodge(0.9),
        show.legend = FALSE,
        na.rm = TRUE # Werte ohne SD (n < 2) erhalten keinen Fehlerbalken
      )
    }
    
    # 8. Ästhetik & ggpubr Integration
    p <- p + scale_fill_manual(values = colourpalette) + 
      scale_color_manual(values = colourpalette)
    
    if (isTRUE(get_input("show_ref_line"))) {
      intercept_val <- if(grepl("Fold|RQ|Relative", yTitle, ignore.case = TRUE)) 1 else 0
      p <- p + geom_hline(yintercept = intercept_val, linetype = "dashed", color = "grey50")
    }
    
    # 9. Signifikanz Labels (Positionierung korrigiert)
    if (isTRUE(get_input("show_significance")) && !is.null(stats_df)) {
      stats_df$signif <- signif_symbol(stats_df$adj.p.value) # same symbols as in the statistics table
      
      df_merged <- merge(df, stats_df[, c("genes", "Label", "signif")], 
                         by.x = c("Gene", "Sample"), by.y = c("genes", "Label"), all.x = TRUE)
      df_merged$signif[is.na(df_merged$signif)] <- "" # e.g. the calibrator/control has no test result

      p <- p + geom_text(
        data = df_merged,
        # ohne SD (n < 2) sitzt das Label direkt über dem Wert
        aes(x = Sample, y = !!sym(PlotDataCol) + dplyr::coalesce(!!sym(PlotDataError), 0), label = signif, group = Gene),
        vjust = -0.5, 
        position = position_dodge(0.9),
        size = axis_text_size * 0.35,
        color = "black",
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
    
    # 10. Theme & Gitterlinien Fix
    p <- p + labs(x = xTitle, y = final_yTitle, title = plotTitle, fill = legendTitle, color = legendTitle) +
      ggpubr::theme_pubr() + 
      theme(
        text = element_text(size = axis_text_size, family = fontfamily),
        plot.title = element_text(size = axis_text_size + 4, face = "bold", hjust = 0.5),
        axis.title = element_text(size = axis_text_size + 2, face = "bold"),
        axis.text.x = element_text(
          angle = xTextAngle, 
          vjust = if(xTextAngle > 0) 0.5 else 1, 
          hjust = if(xTextAngle > 0) 1 else 0.5
        ),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        legend.position = legendPosition,
        legend.title = element_text(face = "bold"),
        strip.text = element_text(size = axis_text_size, face = "bold"),
        panel.grid = element_blank() # Publication ready: keine Linien
      )
    
    if (isTRUE(get_input("facet_by_gene"))) {
      p <- p + facet_wrap(~Gene, scales = "free_y")
    }
    
    if (scalePick != "normal") {
      p <- p + ggpubr::yscale(scalePick, .format = TRUE)
    }
    
    # 11. Ausgabe: ggplot-Objekt für den Export, sonst dasselbe Objekt interaktiv als SVG (ggiraph)
    if (isTRUE(export_mode)) {
      return(p) # Gibt das ggplot-Objekt zurück
    } else {
      return(qrat_girafe(p, width_svg = 10, height_svg = 5))
    }
  }

  # Gemeinsame ggiraph-Einstellungen: Tooltip, Hover-Hervorhebung, Zoom, responsive Breite.
  # Auswahl per Klick nur dort, wo sie ausgewertet wird (Plate View), sonst deaktiviert.
  qrat_girafe <- function(p, width_svg, height_svg, selection = opts_selection(type = "none")) {
    girafe(
      ggobj = p, width_svg = width_svg, height_svg = height_svg,
      options = list(
        opts_sizing(rescale = TRUE),
        opts_tooltip(css = "background:#ffffff;color:#212529;border:1px solid #ced4da;border-radius:6px;padding:6px 10px;font-size:12px;box-shadow:0 2px 6px rgba(0,0,0,0.15);",
                     opacity = 1, use_fill = FALSE),
        opts_hover(css = "opacity:0.7;cursor:pointer;"),
        selection,
        opts_zoom(min = 1, max = 5),
        opts_toolbar(saveaspng = TRUE, pngname = "qRAT_plot")
      )
    )
  }
  
  
  # --- Anwendung im Server-Teil ---
  
  
  # Dieses Reactive bereitet alles für dCq vor (Single Plate)
  prepared_dCq_data <- reactive({
    info <- resultDdct()
    req(info, info$absolute)
    
    # 1. Filtern
    df_plot <- info$absolute %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% input$SamplePicker, Gene %in% input$GenePicker)
    
    validate(need(nrow(df_plot) > 0, "Keine Daten ausgewählt."))
    
    # 2. Spalten-Mapping & Logik
    PlotDataPick <- input$PlotData %||% "RQ"
    col_name_dcq <- if("dCt" %in% names(df_plot)) "dCt" else "dCq"
    col_name_sd  <- if(paste0(col_name_dcq, ".sd") %in% names(df_plot)) paste0(col_name_dcq, ".sd") else "expr.sd"
    
    df_plot[[col_name_dcq]] <- as.numeric(df_plot[[col_name_dcq]])
    df_plot[[col_name_sd]]  <- as.numeric(df_plot[[col_name_sd]])
    if("expr" %in% names(df_plot)) df_plot$expr <- as.numeric(df_plot$expr)
    if("expr.sd" %in% names(df_plot)) df_plot$expr.sd <- as.numeric(df_plot$expr.sd)
    
    if (PlotDataPick == "RQ") {
      col_val <- "expr"; col_sd <- "expr.sd"; ylab <- "Relative Quantity (2^-dCq)"
    } else if (PlotDataPick == "dCq") {
      col_val <- col_name_dcq; col_sd <- col_name_sd; ylab <- "Delta Cq"
    } else {
      df_plot$negdCq <- -(df_plot[[col_name_dcq]])
      col_val <- "negdCq"; col_sd <- col_name_sd; ylab <- "-Delta Cq"
    }
    
    err_col <- if((input$error_type %||% "sd") == "sem") {
      n_val <- if(!is.null(df_plot$n)) as.numeric(df_plot$n) else 3
      df_plot$sem_temp <- df_plot[[col_sd]] / sqrt(n_val)
      "sem_temp"
    } else col_sd
    
    # Wir geben eine Liste zurück, die alles enthält, was der Plotter braucht
    # limma nur rechnen, wenn Signifikanz-Labels angezeigt werden
    return(list(
      plot_df = df_plot,
      stats = if (isTRUE(input$show_significance)) safely(calLimma()) else NULL, # ohne Statistik bleibt der Plot, nur ohne Sterne
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  })
  
  
  output$dcqPlot <- renderGirafe({
    d <- prepared_dCq_data() # Reactive aufrufen
    render_qPCR_plot(d$plot_df, d$stats, d$col_val, d$err_col, d$ylab, input_type = "dCq")
  })
  
  
  
  
  prepared_ddCq_data <- reactive({
    info <- resultDdct()
    req(info, info$relative)
    
    # 1. Filtern basierend auf den DDCt-Pickern
    df_plot <- info$relative %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% input$SamplePickerDDCt, Gene %in% input$GenePickerDDCt)
    
    validate(need(nrow(df_plot) > 0, "Keine Daten ausgewählt."))
    
    # 2. Daten-Typ Wahl
    PlotDataPick <- input$PlotDataDDCt %||% "Fold Change"
    
    # Numerische Konvertierung sicherstellen
    if("FC" %in% names(df_plot)) df_plot$FC <- as.numeric(df_plot$FC)
    if("FC.sd" %in% names(df_plot)) df_plot$FC.sd <- as.numeric(df_plot$FC.sd)
    if("ddCt" %in% names(df_plot)) df_plot$ddCt <- as.numeric(df_plot$ddCt)
    if("ddCt.sd" %in% names(df_plot)) df_plot$ddCt.sd <- as.numeric(df_plot$ddCt.sd)
    
    if (PlotDataPick == "Fold Change") {
      col_val <- "FC"; col_sd <- "FC.sd"; ylab <- "Fold Change (2^-ddCq)"
    } else {
      col_val <- "ddCt"; col_sd <- "ddCt.sd"; ylab <- "Delta Delta Cq"
    }
    
    # 3. Fehlerbalken-Logik (SD vs SEM)
    err_col <- if((input$error_typeDDCt %||% "sd") == "sem") {
      n_val <- if(!is.null(df_plot$n)) as.numeric(df_plot$n) else 3
      df_plot$sem_temp <- df_plot[[col_sd]] / sqrt(n_val)
      "sem_temp"
    } else col_sd
    
    # Rückgabe als Liste für Plot und Export (limma nur bei Signifikanz-Labels)
    return(list(
      plot_df = df_plot,
      stats = if (isTRUE(input$show_significanceDDCt)) safely(calLimma()) else NULL, # ohne Statistik bleibt der Plot, nur ohne Sterne
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  }) 
  
  
  
  
  output$ddctRelGraph2 <- renderGirafe({
    d <- prepared_ddCq_data()
    render_qPCR_plot(d$plot_df, d$stats, d$col_val, d$err_col, d$ylab, input_type = "ddCq")
  })
  
  
  
  
  
  ## ddCt Expression Data Table Single Plate
  
  output$ddctRelative <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdct()
    if (is.null(info)) {
      return(NULL)
    }
    info$relative %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      rename(ddCq = ddCt, ddCq.sd = ddCt.sd)
  })
  
  
  
  
  
  ## MULTIPLE PLATES OUTPUT
  
  
  ## Raw Data Table Multiple Plates
  
  output$multiplePlatesData <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- readMPData()
    if (is.null(info)) NULL else info$data %>% rename(Plate = PlateNumber, Replicate = rp.num)
  })
  
  
  ## Raw Data Table Plate Plot Spatial Original Data Multiple Plates
  
  output$MultiplotCtCard <- renderGirafe({
    info <- readMPData()
    req(info, input$PlateSelect)

    # Zeilennummern beziehen sich auf alle Platten (wie excluded_rows_MP)
    mpv <- info$MPV
    mpv$row_id <- seq_len(nrow(mpv))
    plateData <- mpv %>% filter(PlateNumber == input$PlateSelect)
    req(nrow(plateData) > 0)

    selected <- intersect(excluded_rows_MP(), plateData$row_id)
    plate_view_plot(plateData, row_id = plateData$row_id, selected = as.character(selected))
  })

  # Ausgewählte Wells der angezeigten Platte als manuelle Ausschlüsse übernehmen
  # (Ausschlüsse auf anderen Platten bleiben erhalten)
  observeEvent(input$apply_plate_selection_MP, {
    info <- safely(readMPData())
    req(info, input$PlateSelect)
    plate_rows <- which(info$MPV$PlateNumber == input$PlateSelect)
    selected <- intersect(as.integer(input$MultiplotCtCard_selected), plate_rows)
    excluded_rows_MP(sort(union(setdiff(excluded_rows_MP(), plate_rows), selected)))
    updateRadioButtons(session, "filterModeMP", selected = "manual")
    showNotification(paste0(length(selected), " wells on plate ", input$PlateSelect, " excluded (Manual Filtering, ",
                            length(excluded_rows_MP()), " in total)."), type = "message")
  })

  ## Raw Data Table Cq distribution Original Data Multiple Plates

  output$MultiplotCtDistrib <- renderGirafe({
    info <- multiData()
    req(info)
    cq_density_plot(info$data)
  })
  
  
  
  ## Data Table Multiple Plates Show bad replicates
  
  output$dataMultiplePlatesBadRep <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- multiData()
    req(info)
    data <- info$badReplicates %>%
      select(PlateNumber, Well, Sample, Gene, rp.num, Cq) %>%
      rename(Plate = PlateNumber, Replicate = rp.num)
    if (is.null(info)) NULL else data
  })
  
  
  
  ## Multiple Plates Filter Plot (how many replicates are neglected per sample)
  
  output$MultiplePlatesFilterPlot <- renderGirafe({
    info <- multiData()
    req(info)
    replicate_status_plot(info$badReplicates, info$data)
  })
  
  ## Multiple Plate Boxplot after filtering
  
  output$MultiplePlatesBoxplot <- renderGirafe({
    info <- multiData()
    req(info)
    replicate_boxplot(info$data)
  })
  
  
  
  ## Extracted IPCs Multiple Plates
  
  output$extractedIPCsTable <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    if (input$Id027 == TRUE) {
      info <- extractIPCs()
      data <- info$extrIPCs %>% select(PlateNumber, Well, Sample, Gene, Ct, rp.num) %>%
        rename(Plate = PlateNumber, Cq = Ct, Replicate = rp.num)
      if (is.null(info)) NULL else data
    } else {}
  })
  
  
  
  ## Calibration Factors Table
  
  output$tableCalibrationFactors <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, {
    if (input$Id027 == TRUE) {
      info <- extractIPCs()
      data <- info$CalibrationFactors %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
      if (is.null(info)) NULL else data
    } else {}
  })
  
  ## Inter-Plate Ct values, calibrated with calibration factor
  
  output$interPlateCalibration <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    if (input$Id027 == TRUE) {
      info <- interPlateCalibration()
      data <- info %>%
        select(PlateNumber, Well, Sample, Gene, Ct, rp.num) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
        rename(Plate = PlateNumber, Cq = Ct, Replicate = rp.num)
      if (is.null(info)) NULL else data
    } else {}
  })
  
  
  
  
  
  
  #IPC Comparison Plot
  output$plotIPCcomparison <- renderGirafe({
    req(input$SamplePickerIPCcomparison)
    info <- interPlateCalibration_Plot()
    req(info)
    
    # Daten filtern basierend auf Auswahl
    cal <- info$calibrated %>% 
      filter(Sample_Gene %in% input$SamplePickerIPCcomparison) %>% 
      mutate(Status = "Calibrated")
    
    non <- info$noncalibrated %>% 
      filter(Sample_Gene %in% input$SamplePickerIPCcomparison) %>% 
      mutate(Status = "Raw")
    
    plot_df <- bind_rows(cal, non)
    
    # Faktor-Level festlegen, damit "Raw" immer links und "Calibrated" rechts steht
    plot_df$Status <- factor(plot_df$Status, levels = c("Raw", "Calibrated"))
    
    # Plot erstellen
    p <- ggplot(plot_df, aes(x = Status, y = Ct, color = as.factor(PlateNumber))) +
      # Boxplot für die Streuung
      # (ohne eigene Legende: Raw/Calibrated steht schon auf der x-Achse)
      geom_boxplot(aes(fill = Status), alpha = 0.1, outlier.shape = NA, color = "grey70", show.legend = FALSE) +
      # Einzelne Punkte (Replikate); Hover hebt alle Punkte derselben Platte hervor
      geom_jitter_interactive(width = 0.15, size = 2, alpha = 0.8,
                  aes(tooltip = paste0("<b>", Status, "</b><br>Plate: ", PlateNumber, "<br>Well: ", Well,
                                       "<br>Ct: ", round(Ct, 2), "<br>Gene: ", Gene),
                      data_id = paste("Plate", PlateNumber))) +
      
      # scales = "free_y" erlaubt jedem Gen seine eigene Y-Achse
      facet_wrap(~Sample_Gene, scales = "free_y") +
      theme_minimal() +
      scale_color_brewer(palette = "Set1") + # Bessere Farben für Platten
      labs(title = "Effect of Inter-Plate Calibration",
           x = "", y = "Ct Value", color = "Plate") +
      theme(strip.background = element_rect(fill = "#f8f9fa"),
            legend.position = "bottom")
    
    # Umwandlung in interaktiven Plot
    qrat_girafe(p, width_svg = 9, height_svg = 6)
  })
  
  
  
  
  
  ## dCt Gene Expression Table Multiple Plates
  
  output$ddctAbsoluteMulti <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdctMulti()
    if (is.null(info)) {
      return(NULL)
    }
    info$multiAbsolute %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      rename(dCq = dCt, dCq.sd = dCt.sd, RQ = expr, RQ.sd = expr.sd)
  })
  
  
  
  ## ddCq Gene Expression Table Multiple Plates
  
  output$ddctRelativeMulti <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdctMulti()
    if (is.null(info)) {
      return(NULL)
    }
    info$multiRelative %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      rename(ddCq = ddCt, ddCq.sd = ddCt.sd)
  })
  
  
  
  
  
  
  # # --- Multi Plate Analysis: Absolute dCq Plot ---
  
  
  prepared_MultiAbs_data <- reactive({
    info <- resultDdctMulti()
    req(info, info$multiAbsolute)
    
    # 1. Filtern basierend auf den Multi-Picker-IDs
    df_plot <- info$multiAbsolute %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% input$SamplePickerMulti, Gene %in% input$GenePickerMulti)
    
    validate(need(nrow(df_plot) > 0, "Keine Daten für die Auswahl in Multiple Plate (dCq) vorhanden."))
    
    PlotDataPick <- input$PlotDataMulti %||% "RQ"
    
    # 2. Spalten-Mapping
    col_name_dcq <- if("dCt" %in% names(df_plot)) "dCt" else "dCq"
    col_name_sd  <- if(paste0(col_name_dcq, ".sd") %in% names(df_plot)) paste0(col_name_dcq, ".sd") else "expr.sd"
    
    # Typ-Sicherheit
    df_plot[[col_name_dcq]] <- as.numeric(df_plot[[col_name_dcq]])
    df_plot[[col_name_sd]]  <- as.numeric(df_plot[[col_name_sd]])
    if("expr" %in% names(df_plot)) df_plot$expr <- as.numeric(df_plot$expr)
    if("expr.sd" %in% names(df_plot)) df_plot$expr.sd <- as.numeric(df_plot$expr.sd)
    
    # 3. Auswahl der Werte & Achsenbeschriftung
    if (PlotDataPick == "RQ") {
      col_val <- "expr"; col_sd <- "expr.sd"; ylab <- "Relative Quantity (2^-dCq)"
    } else if (PlotDataPick == "dCq") {
      col_val <- col_name_dcq; col_sd <- col_name_sd; ylab <- "Delta Cq"
    } else {
      df_plot$negdCq <- -(df_plot[[col_name_dcq]])
      col_val <- "negdCq"; col_sd <- col_name_sd; ylab <- "-Delta Cq"
    }
    
    # 4. Error-Berechnung
    err_col <- if((input$error_typeMulti %||% "sd") == "sem") {
      n_val <- if(!is.null(df_plot$n)) as.numeric(df_plot$n) else 3
      df_plot$sem_temp <- df_plot[[col_sd]] / sqrt(n_val)
      "sem_temp"
    } else col_sd
    
    # Rückgabe (limma nur bei Signifikanz-Labels)
    return(list(
      plot_df = df_plot,
      stats = if (isTRUE(input$show_significanceMulti)) safely(calLimmaMulti()) else NULL, # ohne Statistik bleibt der Plot, nur ohne Sterne
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  })
  
  
  output$ddctMultiAbsGraph <- renderGirafe({
    d <- prepared_MultiAbs_data()
    # input_type "MultiAbs" nutzt den Suffix "Multi" für Title, Colors, etc.
    render_qPCR_plot(d$plot_df, d$stats, d$col_val, d$err_col, d$ylab, input_type = "MultiAbs")
  })
  
  
  # --- Multi Plate Analysis: Relative ddCq Plot ---
  
  
  prepared_MultiRel_data <- reactive({
    info <- resultDdctMulti()
    req(info, info$multiRelative)
    
    # 1. Filtern basierend auf den spezifischen Multi-DDCt-Pickern
    df_plot <- info$multiRelative %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% input$SamplePickerDDCtMulti, Gene %in% input$GenePickerDDCtMulti)
    
    validate(need(nrow(df_plot) > 0, "Keine Daten für die Auswahl in Multiple Plate (ddCq) vorhanden."))
    
    # 2. Daten-Wahl
    PlotDataPick <- input$PlotDataDDCtMulti %||% "Fold Change"
    
    # Typ-Sicherheit
    if("FC" %in% names(df_plot)) df_plot$FC <- as.numeric(df_plot$FC)
    if("FC.sd" %in% names(df_plot)) df_plot$FC.sd <- as.numeric(df_plot$FC.sd)
    if("ddCt" %in% names(df_plot)) df_plot$ddCt <- as.numeric(df_plot$ddCt)
    if("ddCt.sd" %in% names(df_plot)) df_plot$ddCt.sd <- as.numeric(df_plot$ddCt.sd)
    
    if (PlotDataPick == "Fold Change") {
      col_val <- "FC"; col_sd <- "FC.sd"; ylab <- "Fold Change (2^-ddCq)"
    } else {
      col_val <- "ddCt"; col_sd <- "ddCt.sd"; ylab <- "Delta Delta Cq"
    }
    
    # 3. Error-Berechnung
    err_col <- if((input$error_typeDDCtMulti %||% "sd") == "sem") {
      n_val <- if(!is.null(df_plot$n)) as.numeric(df_plot$n) else 3
      df_plot$sem_temp <- df_plot[[col_sd]] / sqrt(n_val)
      "sem_temp"
    } else col_sd

    # limma nur bei Signifikanz-Labels
    return(list(
      plot_df = df_plot,
      stats = if (isTRUE(input$show_significanceDDCtMulti)) safely(calLimmaMulti()) else NULL, # ohne Statistik bleibt der Plot, nur ohne Sterne
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  })  
  
  
  output$ddctRelGraphMulti <- renderGirafe({
    d <- prepared_MultiRel_data()
    # Nutzt input_type "MultiRel" -> Suffix "DDCtMulti"
    render_qPCR_plot(d$plot_df, d$stats, d$col_val, d$err_col, d$ylab, input_type = "MultiRel")
  })  
  
  
  
  ## MIQE Checks
  ## Single Plate MIQE check
  
  # Sample/Gene-Kombinationen mit < 2 gültigen technischen Replikaten (Tabelle + Statusbox teilen sich das Reactive)
  lowReplicatesSP <- reactive({
    info <- setData()
    req(info)
    info$data %>%
      filter(!is.na(Ct)) %>%
      count(Sample, Gene) %>%
      filter(n < 2)
  })

  output$countTable <- renderDT(lowReplicatesSP())

  output$SinglePlateMIQEcheck <- renderUI({
    countReplicates <- lowReplicatesSP()

    if (dim(countReplicates)[1] == 0) {
      tagList(
        withTags({
          # style(".fa-check {color:#E87722}")
          div(
            class = "bs-component",
            div(
              class = "alert alert-success",
              div(class = "float-right", icon("check", "fa-3x")),
              h4(class = "alert-heading", "Number of technical replicates >=2")
            )
          )
        }),
      )
    } else {
      tagList(
        withTags({
          div(
            class = "bs-component",
            div(
              class = "alert alert-warning",
              div(class = "float-right", icon("xmark", "fa-3x")),
              h4(class = "alert-heading", "Replicates")
            ),
            p(div(class = "border border-warning", strong("Following samples are below recommended number of technical replicates:"), DTOutput("countTable")))
          )
        })
      )
    }
  })
  
  ## Multiple Plate MIQE check
  
  lowReplicatesMP <- reactive({
    info <- multiData()
    req(info)
    info$data %>%
      filter(!is.na(Ct)) %>%
      count(Sample, Gene) %>%
      filter(n < 2)
  })

  # eigene Output-ID, damit Single- und Multi-Plate sich die Tabelle nicht gegenseitig überschreiben
  output$countTableMulti <- renderDT(lowReplicatesMP())

  output$MultiplePlateMIQEcheck <- renderUI({
    countReplicates <- lowReplicatesMP()

    if (dim(countReplicates)[1] == 0) {
      tagList(
        withTags({
          # style(".fa-check {color:#E87722}")
          div(
            class = "bs-component",
            div(
              class = "alert alert-success",
              div(class = "float-right", icon("check", "fa-3x")),
              h4(class = "alert-heading", "Number of technical replicates >=2")
            )
          )
        }),
      )
    } else {
      tagList(
        withTags({
          div(
            class = "bs-component",
            div(
              class = "alert alert-warning",
              div(class = "float-right", icon("xmark", "fa-3x")),
              h4(class = "alert-heading", "Replicates")
            ),
            p(div(class = "border border-warning", strong("Following samples are below recommended number of technical replicates:"), DTOutput("countTableMulti")))
          )
        })
      )
    }
  })



  ## Limma Outputs
  
  output$resultLimma <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    dt <- calLimma()
    # Hinweis statt leerer Tabelle (z. B. nur eine Probe, kein Vergleich gewählt, keine Referenzgene)
    validate(need(!is.null(dt), "No statistics available: select reference gene(s) and a comparison between at least two samples."))
    
    dt %>%
      select(Label, genes, t.test, p.value, adj.p.value) %>%
      rename(Sample = Label, Gene = genes) %>%
      mutate(significance = signif_symbol(adj.p.value)) %>% # same symbols as in the plots (qPCR.R)
      mutate(across(c("p.value", "adj.p.value"), label_pvalue(accuracy = 0.0001))) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  })
  
  output$resultLimmaMulti <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    dt <- calLimmaMulti()
    # Hinweis statt leerer Tabelle (z. B. nur eine Probe, kein Vergleich gewählt, keine Referenzgene)
    validate(need(!is.null(dt), "No statistics available: select reference gene(s) and a comparison between at least two samples."))
    dt %>%
      select(Label, genes, t.test, p.value, adj.p.value) %>%
      rename(Sample = Label, Gene = genes) %>%
      mutate(significance = signif_symbol(adj.p.value)) %>% # same symbols as in the plots (qPCR.R)
      mutate(across(c("p.value", "adj.p.value"), label_pvalue(accuracy = 0.0001))) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  })
  
  
  
  
  ## conditions for conditionalpanel Single Plate & Multiple Plate Upload
  output$fileUploadedSingle <- reactive({
    return(!is.null(readSinglePlateData()))
  })
  output$fileUploadedMulti <- reactive({
    return(!is.null(multiData()))
  })
  
  lapply(c("fileUploadedMulti", "fileUploadedSingle"), function(x) outputOptions(output, x, suspendWhenHidden = FALSE))
  
  
  
  
  ### Comparison Boxes for Statistical Analysis
  ### Single Plate
  # Track the number of input boxes to render
  counter <- reactiveValues(n = 0)
  
  # Track all user inputs
  AllInputs <- reactive({
    x <- reactiveValuesToList(input)
  })
  
  observeEvent(input$add_btn, {
    counter$n <- counter$n + 1
  })
  observeEvent(input$rm_btn, {
    if (counter$n > 0) counter$n <- counter$n - 1
  })
  
  comparisonBoxes <- reactive({
    n <- counter$n
    info <- setData()
    
    if (n > 0) {
      isolate({
        lapply(seq_len(n), function(i) {
          virtualSelectInput(
            inputId = paste0("CompsMultiple", i),
            label = paste0("Comparison ", i),
            choices = info$CompSamples,
            multiple = TRUE,
            maxValues = 2,
            disableAllOptionsSelectedText = TRUE,
            showValueAsTags = TRUE,
            selected = AllInputs()[[paste0("CompsMultiple", i)]]
          )
        })
      })
    }
  })
  
  # Paarweise Vergleiche aus den Vergleichsboxen ("A:B,C:D", Leerzeichen -> ".").
  # Hängt direkt von der Auswahl in den Boxen ab; vorher stand sie in isolate(), sodass eine
  # geänderte Auswahl limma nicht neu rechnen ließ. Nur vollständige Paare (2 Proben) zählen,
  # damit eine halb ausgefüllte Box keine fehlerhafte Kontrast-Formel erzeugt.
  build_paired_comparisons <- function(n, prefix) {
    if (n == 0) return(NULL)
    pairs <- lapply(seq_len(n), function(i) input[[paste0(prefix, i)]])
    pairs <- Filter(function(x) length(x) == 2, pairs)
    if (length(pairs) == 0) return(NULL)
    multiComps <- paste(vapply(pairs, paste, character(1), collapse = ":"), collapse = ",")
    gsub(" ", ".", multiComps)
  }

  listOfCompsMultiple <- reactive({
    build_paired_comparisons(counter$n, "CompsMultiple")
  })
  
  output$comparisonBoxes_ui <- renderUI({
    comparisonBoxes()
  })
  
  ### Comparison Boxes for Statistical Analysis
  ### Multiple Plates
  # Track the number of input boxes to render
  counterM <- reactiveValues(n = 0)

  # (AllInputs ist oben bei Single Plate definiert)

  observeEvent(input$add_btnM, {
    counterM$n <- counterM$n + 1
  })
  observeEvent(input$rm_btnM, {
    if (counterM$n > 0) counterM$n <- counterM$n - 1
  })
  
  comparisonBoxesM <- reactive({
    n <- counterM$n
    info <- multiData2()
    
    if (n > 0) {
      isolate({
        lapply(seq_len(n), function(i) {
          virtualSelectInput(
            inputId = paste0("CompsMultipleM", i),
            label = paste0("Comparison ", i),
            choices = info$CompSamples,
            multiple = TRUE,
            maxValues = 2,
            disableAllOptionsSelectedText = TRUE,
            showValueAsTags = TRUE,
            selected = AllInputs()[[paste0("CompsMultipleM", i)]]
          )
        })
      })
    }
  })
  
  listOfCompsMultipleM <- reactive({
    build_paired_comparisons(counterM$n, "CompsMultipleM")
  })
  
  output$comparisonBoxesM_ui <- renderUI({
    comparisonBoxesM()
  })
  
  
  
  # --- 5. NAVIGATION & UI-HELPER ---
  
  observeEvent(input$link_to_tab_single, {
    updateNavbarPage(session, "tabs", selected = "Single Plate Analysis")
  })
  
  observeEvent(input$link_to_tab_multi, {
    updateNavbarPage(session, "tabs", selected = "Multiple Plate Analysis")
  })
  
  
  # Zuerst den Output definieren
  output$fileUploaded <- reactive({
    !is.null(input$dtfile)
  })
  
  # Dann die Optionen setzen
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  
  
  
  
  
  # --- MANUAL FILTERING LOGIC Single Plate---
  
  # Speichert die IDs (oder Zeilennummern), die manuell ausgeschlossen wurden
  excluded_rows <- reactiveVal(numeric(0))
  
  # Reset, wenn eine neue Datei geladen wird
  observeEvent(input$dtfile, {
    excluded_rows(numeric(0))
  })
  
  # 1. Modal öffnen
  observeEvent(input$open_manual_filter, {
    req(input$dtfile)
    info <- safely(readSinglePlateData())
    df <- info$data
    if(is.null(df)) return()
    
    # Füge eine temporäre ID hinzu, um Zeilen eindeutig zu identifizieren
    df$RowID <- 1:nrow(df)
    
    showModal(modalDialog(
      title = "Manual Outlier Selection",
      size = "l", # Large modal
      h5("Click on rows to mark them as outliers (excluded)."),
      DTOutput("manual_filter_table"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_manual_filter", "Apply Exclusions", class = "btn-success")
      )
    ))
  })
  
  # 2. Tabelle im Modal rendern
  output$manual_filter_table <- renderDT({
    info <- readSinglePlateData()
    df <- info$data
    req(df)
    
    # Wir zeigen Sample, Gene, Ct und Position (falls vorhanden) an
    datatable(df %>% select(any_of(c("Sample", "Gene", "Ct", "Cq", "Pos", "Well"))),
              selection = list(mode = 'multiple', selected = excluded_rows()),
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 3. Speichern der Auswahl
  observeEvent(input$save_manual_filter, {
    # Die Indices der in der Tabelle ausgewählten Zeilen holen
    selected_indices <- input$manual_filter_table_rows_selected
    
    # Speichern
    excluded_rows(selected_indices)
    
    removeModal()
    
    # Optional: Eine Notification anzeigen
    showNotification(paste(length(selected_indices), "data points manually excluded."), type = "message")
  })
  
  # UI Feedback: Wie viele sind ausgeschlossen?
  output$manual_filter_status <- renderUI({
    count <- length(excluded_rows())
    if(count > 0) {
      tags$div(style="color: red; font-weight: bold;", 
               paste(count, "data points currently excluded."))
    } else {
      tags$div(style="color: green;", "No manual exclusions.")
    }
  })  
  
  
  
  
  # --- MANUAL FILTERING LOGIC (MULTIPLE PLATES) ---
  
  # Speicher für ausgeschlossene Zeilen (MP)
  excluded_rows_MP <- reactiveVal(numeric(0))
  
  # Reset, wenn neue Platten geladen werden
  observeEvent(input$plates, {
    excluded_rows_MP(numeric(0))
  })
  
  # 1. Modal öffnen (MP)
  observeEvent(input$open_manual_filter_MP, {
    req(input$plates)
    info <- safely(readMPData())
    df <- info$data
    if(is.null(df)) return()
    
    showModal(modalDialog(
      title = "Manual Outlier Selection (Multiple Plates)",
      size = "l", 
      h5("Click on rows to mark them as outliers (excluded)."),
      DTOutput("manual_filter_table_MP"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_manual_filter_MP", "Apply Exclusions", class = "btn-success")
      )
    ))
  })
  
  # 2. Tabelle im Modal rendern (MP)
  output$manual_filter_table_MP <- renderDT({
    info <- readMPData()
    df <- info$data
    req(df)
    
    # Wir zeigen PlateNumber, Sample, Gene, Ct an
    cols_to_show <- c("PlateNumber", "Sample", "Gene", "Ct", "Cq", "Pos", "Well")
    
    datatable(df %>% select(any_of(cols_to_show)),
              selection = list(mode = 'multiple', selected = excluded_rows_MP()),
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 3. Speichern der Auswahl (MP)
  observeEvent(input$save_manual_filter_MP, {
    selected_indices <- input$manual_filter_table_MP_rows_selected
    excluded_rows_MP(selected_indices)
    removeModal()
    showNotification(paste(length(selected_indices), "data points manually excluded (MP)."), type = "message")
  })
  
  # Status Anzeige UI (MP)
  output$manual_filter_status_MP <- renderUI({
    count <- length(excluded_rows_MP())
    if(count > 0) {
      tags$div(style="color: red; font-weight: bold;", 
               paste(count, "data points currently excluded."))
    } else {
      tags$div(style="color: green;", "No manual exclusions.")
    }
  })    
  
  
  
  
  
  # --- MANUELLE PROBEN-GRUPPIERUNG SP ---
  manual_group_cnt <- reactiveVal(0)
  active_group_ids <- reactiveVal(c())
  confirmedGrouping <- reactiveVal(NULL)
  
  createGroupUI <- function(id, name = "", selected_samples = NULL, available_samples = c("")) {
    tags$div(
      id = paste0("group_container_", id),
      style = "border: 1px solid #dee2e6; padding: 15px; margin-bottom: 10px; border-radius: 8px; background: #fdfdfd;",
      layout_column_wrap(
        width = 1/2,
        textInput(paste0("man_group_name_", id), "Group Name", value = name, placeholder = "e.g. Treatment-A"),
        virtualSelectInput(paste0("man_group_samples_", id), "Assign Samples", 
                           choices = available_samples, selected = selected_samples,
                           multiple = TRUE, search = TRUE)
      ),
      tags$div(style = "text-align: right;",
               actionButton(paste0("remove_grp_btn_", id), "Remove Group", 
                            icon = icon("trash-can"), class = "btn-link text-danger btn-sm"))
    )
  }
  
  observeEvent(input$open_grouping_modal, {
    showModal(modalDialog(
      title = "Configure Biological Replicates (Single Plate)",
      size = "l", easyClose = FALSE,
      footer = tagList(modalButton("Cancel"), actionButton("save_grouping_btn", "Apply & Save", class = "btn-success")),
      tags$div(id = "manual_groups_wrapper", style = "min-height: 100px;"),
      hr(),
      actionButton("add_manual_group_btn", "Add New Group", icon = icon("plus-circle"), class = "btn-outline-primary")
    ))
    
    current_data <- confirmedGrouping()
    available_samples <- c("")
    try({
      raw_info <- readSinglePlateData()
      if (!is.null(raw_info)) available_samples <- sort(unique(as.character(raw_info$data$Sample)))
    }, silent = TRUE)
    
    if (!is.null(current_data) && length(current_data) > 0) {
      saved_ids <- names(current_data)
      active_group_ids(saved_ids)
      for (id in saved_ids) {
        insertUI(selector = "#manual_groups_wrapper", where = "beforeEnd",
                 ui = createGroupUI(id, current_data[[id]]$name, current_data[[id]]$samples, available_samples))
      }
    }
  })
  
  observeEvent(input$add_manual_group_btn, {
    new_id <- as.character(manual_group_cnt() + 1)
    manual_group_cnt(as.numeric(new_id))
    active_group_ids(c(active_group_ids(), new_id))
    
    available_samples <- c("")
    try({
      raw_info <- readSinglePlateData()
      if (!is.null(raw_info)) available_samples <- sort(unique(as.character(raw_info$data$Sample)))
    }, silent = TRUE)
    
    insertUI(selector = "#manual_groups_wrapper", where = "beforeEnd",
             ui = createGroupUI(new_id, available_samples = available_samples))
  })
  
  observeEvent(input$save_grouping_btn, {
    group_data <- list()
    for (id in active_group_ids()) {
      name <- input[[paste0("man_group_name_", id)]]
      samples <- input[[paste0("man_group_samples_", id)]]
      if (!is.null(name) && name != "" && length(samples) > 0) {
        group_data[[id]] <- list(name = name, samples = samples)
      }
    }
    confirmedGrouping(group_data)
    removeModal()
    showNotification("Grouping updated.", type = "message")
  })
  
  # Genau ein Remove-Observer pro Gruppen-ID. Vorher entstanden bei jeder Änderung neue Kopien,
  # und ohne local() hörten alle Observer nur auf den Button der zuletzt angelegten Gruppe.
  observed_group_ids <- character(0)
  observe({
    new_ids <- setdiff(active_group_ids(), observed_group_ids)
    observed_group_ids <<- c(observed_group_ids, new_ids)
    for (new_id in new_ids) {
      local({
        id <- new_id
        observeEvent(input[[paste0("remove_grp_btn_", id)]], {
          removeUI(selector = paste0("#group_container_", id))
          active_group_ids(active_group_ids()[active_group_ids() != id])
        }, ignoreInit = TRUE)
      })
    }
  })
  
  
  # --- MANUELLE PROBEN-GRUPPIERUNG MP ---
  manual_group_cnt_multi <- reactiveVal(0)
  active_group_ids_multi <- reactiveVal(c())
  confirmedGroupingMulti <- reactiveVal(NULL)
  
  # Hilfsfunktion (IDs sind spezifisch für Multi-Plate)
  createGroupUIMulti <- function(id, name = "", selected_samples = NULL, available_samples = c("")) {
    tags$div(
      id = paste0("group_container_multi_", id),
      style = "border: 1px solid #dee2e6; padding: 15px; margin-bottom: 10px; border-radius: 8px; background: #fdfdfd;",
      layout_column_wrap(
        width = 1/2,
        textInput(paste0("man_group_name_multi_", id), "Group Name", value = name, placeholder = "e.g. Wildtype-Group"),
        virtualSelectInput(paste0("man_group_samples_multi_", id), "Assign Samples", 
                           choices = available_samples, selected = selected_samples,
                           multiple = TRUE, search = TRUE)
      ),
      tags$div(style = "text-align: right;",
               actionButton(paste0("remove_grp_btn_multi_", id), "Remove Group", 
                            icon = icon("trash-can"), class = "btn-link text-danger btn-sm"))
    )
  }
  
  # 1. Modal öffnen
  observeEvent(input$open_grouping_modalMulti, {
    
    mp_info <- safely(multiData())
    req(mp_info) 
    
    removeModal() 
    
    showModal(modalDialog(
      title = "Configure Biological Replicates (Multiple Plates)",
      size = "l", 
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"), 
        actionButton("save_grouping_btn_multi", "Apply & Save", class = "btn-success")
      ),
      tags$div(id = "manual_groups_wrapper_multi", style = "min-height: 100px;"),
      hr(),
      actionButton("add_manual_group_btn_multi", "Add New Group", icon = icon("plus-circle"), class = "btn-outline-primary")
    ))
    
    # Probenliste aus der multiData Reactive beziehen
    available_samples <- c("")
    if (!is.null(mp_info)) {
      # Wir greifen auf die Spalte 'Sample' im Daten-Element der Liste zu
      available_samples <- sort(unique(as.character(mp_info$data$Sample)))
    }
    
    # Vorhandene Gruppen wiederherstellen
    current_data <- confirmedGroupingMulti()
    if (!is.null(current_data) && length(current_data) > 0) {
      saved_ids <- names(current_data)
      active_group_ids_multi(saved_ids)
      for (id in saved_ids) {
        insertUI(selector = "#manual_groups_wrapper_multi", where = "beforeEnd",
                 ui = createGroupUIMulti(id, current_data[[id]]$name, current_data[[id]]$samples, available_samples))
      }
    }
  })
  
  # 2. Neue Gruppe hinzufügen
  observeEvent(input$add_manual_group_btn_multi, {
    new_id <- as.character(manual_group_cnt_multi() + 1)
    manual_group_cnt_multi(as.numeric(new_id))
    active_group_ids_multi(c(active_group_ids_multi(), new_id))
    
    available_samples <- c("")
    mp_info <- safely(multiData())
    if (!is.null(mp_info)) {
      available_samples <- sort(unique(as.character(mp_info$data$Sample)))
    }
    
    insertUI(selector = "#manual_groups_wrapper_multi", where = "beforeEnd",
             ui = createGroupUIMulti(new_id, available_samples = available_samples))
  })
  
  # 3. Speichern
  observeEvent(input$save_grouping_btn_multi, {
    group_data <- list()
    for (id in active_group_ids_multi()) {
      name <- input[[paste0("man_group_name_multi_", id)]]
      samples <- input[[paste0("man_group_samples_multi_", id)]]
      if (!is.null(name) && name != "" && length(samples) > 0) {
        group_data[[id]] <- list(name = name, samples = samples)
      }
    }
    confirmedGroupingMulti(group_data)
    removeModal()
    showNotification("Multi-Plate Grouping updated.", type = "message")
  })
  
  # 4. Löschen von Gruppen (genau ein Observer pro Gruppen-ID, siehe Single Plate)
  observed_group_ids_multi <- character(0)
  observe({
    new_ids <- setdiff(active_group_ids_multi(), observed_group_ids_multi)
    observed_group_ids_multi <<- c(observed_group_ids_multi, new_ids)
    for (new_id in new_ids) {
      local({
        id <- new_id
        observeEvent(input[[paste0("remove_grp_btn_multi_", id)]], {
          removeUI(selector = paste0("#group_container_multi_", id))
          active_group_ids_multi(active_group_ids_multi()[active_group_ids_multi() != id])
        }, ignoreInit = TRUE)
      })
    }
  })
  
  
  
  #Reference Gene Finder
  # 1. Dynamische Auswahl der Reference-Genes Single Plate
  observe({
    # Wir nehmen die Daten nach dem Filtering
    info <- safely(setData())
    req(info)
    available_genes <- sort(unique(as.character(info$data$Gene)))
    
    updateVirtualSelect("ref_candidates", choices = available_genes)
  })
  
  # Dynamische Auswahl der Reference-Genes Multi Plate
  observe({
    info <- safely(multiData2())
    req(info)
    # Wir nehmen die Gene aus dem bereinigten Multi-Plate Datensatz
    available_genes <- sort(unique(as.character(info$SamplesNoIPCs$Gene)))
    updateVirtualSelect("ref_candidates_multi", choices = available_genes)
  })
  
  
  # 2. Berechnung Triggern Single Plate
  validationResults <- eventReactive(input$run_val_btn, {
    req(input$ref_candidates)
    
    # 1. Daten holen & Filtern
    plot_data <- setData()$data %>%
      filter(Gene %in% input$ref_candidates) %>%
      group_by(Sample, Gene) %>%
      summarise(Ct = mean(Ct, na.rm = TRUE), .groups = "drop")
    
    # 2. In Matrix-Format bringen
    wide_df <- plot_data %>%
      pivot_wider(names_from = Gene, values_from = Ct)
    
    # Umwandlung in Matrix: Erste Spalte (Sample) als Row-Names nutzen
    ct_mat <- as.matrix(wide_df[, -1]) # Alles außer der ersten Spalte
    rownames(ct_mat) <- wide_df$Sample  # Die erste Spalte als Rownames setzen
    
    # 3. Algorithmus ausführen
    withProgress(message = 'Calculating Stability...', value = 0.5, {
      res <- run_genorm_logic(ct_mat)
      setProgress(1)
    })
    
    return(res)
  })
  
  # Plot für die Stabilität
  output$ref_stability_plot <- renderGirafe({
    res <- validationResults()
    req(res)
    
    p <- ggplot(res, aes(x = reorder(Gene, M), y = M, fill = M,
                         tooltip = paste0("<b>", Gene, "</b><br>M value: ", round(M, 4), "<br>Rank: ", Rank),
                         data_id = Gene)) +
      geom_bar_interactive(stat = "identity", alpha = 0.8) +
      scale_fill_gradient(low = "#2ecc71", high = "#e74c3c") +
      theme_minimal() +
      labs(x = "Candidate Genes", y = "geNorm M value (lower = more stable)",
           title = "geNorm Ranking") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    qrat_girafe(p, width_svg = 7, height_svg = 5)
  })
  
  # 3. Output: Ranking Tabelle
  output$ref_ranking_table <- renderDT({
    res <- validationResults()
    req(res)
    
    datatable(res, 
              extensions = c("Buttons"),
              options = list(
                dom = 'tB', # 't' für Tabelle, 'B' für Buttons
                buttons = c('copy', 'csv', 'excel'),
                pageLength = 15
              ), 
              rownames = FALSE) %>%
      formatRound("M", 3)
  }) 
  
  
  ##Reference Gene Finder Multiplate
  validationResultsMulti <- eventReactive(input$run_val_btn_multi, {
    req(multiData2())
    candidates <- input$ref_candidates_multi
    req(candidates)
    
    # bei aktiver IPC-Kalibrierung die kalibrierten Cq-Werte (wie in ddCq/limma), sonst die unkalibrierten
    raw_df <- multiData2()$SamplesCalibrated
    num_plates <- length(unique(raw_df$PlateNumber))
    
    # --- PRÜFUNG: Gene auf allen Platten vorhanden ---
    gene_presence <- raw_df %>%
      filter(Gene %in% candidates) %>%
      group_by(Gene) %>%
      summarise(n_plates = n_distinct(PlateNumber))
    
    incomplete_genes <- gene_presence %>% 
      filter(n_plates < num_plates) %>% 
      pull(Gene)
    
    if (length(incomplete_genes) > 0) {
      showModal(modalDialog(
        title = "Warning: Incomplete Gene Data",
        tags$p("The following genes are not present on all plates:"),
        tags$ul(lapply(incomplete_genes, tags$li)),
        tags$p("This can lead to biased stability values. We recommend using genes that were measured across all plates."),
        footer = modalButton("I understand")
      ))
    }
    
    # --- EIGENTLICHE BERECHNUNG ---
    plot_data <- raw_df %>%
      filter(Gene %in% candidates) %>%
      group_by(Sample, Gene) %>%
      summarise(Ct = mean(Ct, na.rm = TRUE), .groups = "drop")
    
    # Pivot & Matrix
    wide_df <- plot_data %>%
      pivot_wider(names_from = Gene, values_from = Ct)
    
    ct_mat <- as.matrix(wide_df[, -1])
    rownames(ct_mat) <- wide_df$Sample
    
    withProgress(message = 'Calculating Multi-Plate Stability...', value = 0.5, {
      res <- run_genorm_logic(ct_mat)
      setProgress(1)
    })

    # Abdeckung hinzufügen
    coverage <- gene_presence %>% 
      mutate(Coverage = paste0(round((n_plates / num_plates) * 100), "%")) %>%
      select(Gene, Coverage)
    
    res <- res %>% left_join(coverage, by = "Gene")
    
    return(res)
  })
  
  
  
  output$ref_stability_plot_multi <- renderGirafe({
    res <- validationResultsMulti()
    req(res)
    
    # Wir erstellen den Plot analog zum Single-Plate Plot, fügen aber die Coverage-Information hinzu.
    p <- ggplot(res, aes(x = reorder(Gene, M), y = M, fill = M,
                         tooltip = paste0("<b>", Gene, "</b>",
                                          "<br>M value: ", round(M, 4),
                                          "<br>Plate Coverage: ", Coverage),
                         data_id = Gene)) +
      geom_bar_interactive(stat = "identity", alpha = 0.8) +
      # Farbskala: Grün (stabil) zu Rot (instabil)
      scale_fill_gradient(low = "#2ecc71", high = "#e74c3c") +
      theme_minimal() +
      labs(x = "Candidate Genes", 
           y = "geNorm M value (lower = more stable)",
           title = "Multi-Plate geNorm Ranking") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    qrat_girafe(p, width_svg = 7, height_svg = 5)
  })
  
  output$ref_ranking_table_multi <- renderDT({
    res <- validationResultsMulti()
    req(res)
    
    datatable(res, rownames = FALSE) %>%
      formatRound("M", 3) %>%
      formatStyle(
        'Coverage',
        color = styleEqual("100%", "green", "red"),
        fontWeight = styleEqual("100%", "normal", "bold")
      )
  })
  
  
  
  
  # --- NORMFINDER LOGIK ---
  
  # --- geNorm (Vandesompele et al. 2002, umgesetzt in ctrlGene) ---
  # ct_matrix: Zeilen = Proben, Spalten = Kandidatengene (mittlere Cq je Probe).
  # Rückgabe: Gene, M (geNorm-Stabilitätswert, kleiner = stabiler), Rank (das stabilste Paar teilt Rang 1)
  run_genorm_logic <- function(ct_matrix) {
    # geNorm vergleicht Gene paarweise über alle Proben und braucht daher vollständige Proben
    ct_matrix <- ct_matrix[stats::complete.cases(ct_matrix), , drop = FALSE]
    validate(
      need(ncol(ct_matrix) >= 2, "geNorm needs at least two candidate genes."),
      need(nrow(ct_matrix) >= 2, "geNorm needs at least two samples with values for all candidate genes.")
    )
    res <- ctrlGene::geNorm(ct_matrix, ctVal = TRUE)
    # ctrlGene listet die schrittweise ausgeschlossenen Gene (instabilstes zuerst) und zuletzt das stabilste Paar
    # ("A-B"); das Paar wird über die Spaltennamen bestimmt, damit Gennamen mit "-" kein Problem sind
    singles <- as.character(res$Genes[-nrow(res)])
    pair <- setdiff(colnames(ct_matrix), singles)
    m_pair <- res$Avg.M[nrow(res)]
    m_singles <- rev(res$Avg.M[-nrow(res)])
    data.frame(
      Gene = c(pair, rev(singles)),
      M = c(rep(m_pair, length(pair)), m_singles),
      Rank = c(rep(1L, length(pair)), seq_along(singles) + 1L)
    )
  }
  
  
  #Reference Gene Finder Apply Button Single Plate
  observeEvent(input$apply_best_ref, {
    res <- safely(validationResults())
    req(res)
    
    # Nimm die Top 2 Gene
    best_genes <- res$Gene[1:2]
    
    # Update das Haupt-Input für Referenzgene (aus dem Analyse-Tab)
    updateVirtualSelect("Refs", selected = best_genes)
    
    showNotification(paste("Applied", paste(best_genes, collapse = ", "), "as reference genes."), type = "message")
  }) 
  
  #Reference Gene Finder Apply Button Multi Plate
  observeEvent(input$apply_best_ref_multi, {
    res <- safely(validationResultsMulti())
    req(res)
    
    # Die Top 2 stabilsten Gene
    best_genes <- res$Gene[1:2]
    
    # Update das VirtualSelect im Multi-Plate Tab
    updateVirtualSelect("RefsM", selected = best_genes)
    
    showNotification(
      paste("Multi-Plate: Applied", paste(best_genes, collapse = ", "), "as reference genes."), 
      type = "message"
    )
  })
  
  
  
  # Erklärung des Reference Finders (Single und Multi Plate)
  genorm_description <- function() {
    tagList(
      tags$h5("geNorm Algorithm", style = "color: #2c3e50;"),
      tags$p("geNorm ranks candidate reference genes by their expression stability across all samples."),

      tags$div(
        style = "background-color: #f8f9fa; border-left: 5px solid #2ecc71; padding: 15px; margin: 10px 0;",
        tags$b("Interpretation:"),
        tags$ul(
          tags$li(tags$b("M value:"), " average pairwise variation of a gene with all other candidate genes."),
          tags$li(tags$b("Rule of Thumb:"), " A lower M value indicates higher expression stability."),
          tags$li("Values ", tags$b("< 0.5"), " are generally considered stable in homogeneous samples (up to 1 in heterogeneous samples).")
        )
      ),

      tags$h6("How it works:"),
      tags$p("For every candidate, geNorm computes the standard deviation of its log-ratios with each other candidate across the samples; the mean of these values is the M value. The least stable gene is removed step by step and M is recalculated, until the two most stable genes remain (shared rank 1). The calculation uses the mean Cq per sample and gene; samples without a value for all candidate genes are skipped."),

      tags$hr(),

      tags$div(
        style = "font-size: 0.85rem; color: #7f8c8d;",
        tags$p(tags$i("Reference: Vandesompele J, De Preter K, Pattyn F, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol. 2002;3(7):research0034."))
      )
    )
  }

  output$algo_description <- renderUI(genorm_description())

  output$algo_description_multi <- renderUI(genorm_description())
  
  
  
  
  
  #Info-Box zu IPC
  output$ipc_info_text <- renderUI({
    tagList(
      tags$h6("The Ct-Shift Principle"),
      tags$p("Since Ct values are logarithmic (log2), technical differences between plates manifest as a ", 
             tags$b("constant numerical offset"), " rather than a ratio."),
      
      tags$div(
        style = "background-color: #e9f7ef; border-left: 4px solid #2ecc71; padding: 10px; margin-bottom: 10px;",
        tags$p(style = "margin-bottom: 5px;", tags$b("Calculation:")),
        tags$code("Ct(corrected) = Ct(raw) - (Plate_IPC_Mean - Global_IPC_Mean)")
      ),
      
      tags$p("By subtracting the plate-specific deviation (Mean), we align the baseline of all plates without distorting the biological fold-change between samples."),
      
      tags$h6("Calibration Check Plot"),
      tags$p("The plot shows the distribution of Ct values for selected samples across all plates. In an ideal calibration, the 'Calibrated' points should cluster much tighter than the 'Raw' points.")
    )
  })
  
  
  
  
  #Plot Export Download Handler
  create_qPCR_downloadHandler <- function(input, data_reactive, type_suffix, col_name, err_name, y_lab) {
    downloadHandler(
      filename = function() {
        # Nutzt das Format-Input aus der Sidebar (z.B. input$plot_format_dCq)
        format <- input[[paste0("plot_format_", type_suffix)]] %||% "png"
        paste0("qPCR_", type_suffix, "_Plot_", Sys.Date(), ".", format)
      },
      content = function(file) {
        req(data_reactive())
        d <- data_reactive() # Das ist jetzt unsere Liste von oben
        
        p_static <- render_qPCR_plot(
          df = d$plot_df,
          stats_df = d$stats,
          PlotDataCol = d$col_val, # Nimmt dynamisch dCq oder RQ
          PlotDataError = d$err_col,
          yTitle = d$ylab,
          input_type = type_suffix,
          export_mode = TRUE 
        )
        
        # Export-Parameter dynamisch abgreifen
        w   <- input[[paste0("plot_w_", type_suffix)]] %||% 7
        h   <- input[[paste0("plot_h_", type_suffix)]] %||% 5
        res_scale <- input[[paste0("plot_res_", type_suffix)]] %||% 3
        fmt <- input[[paste0("plot_format_", type_suffix)]] %||% "png"
        
        ggsave(
          filename = file,
          plot = p_static,
          device = fmt,
          width = w,
          height = h,
          units = "in",
          dpi = 100 * res_scale
        )
      }
    )
  }
  
  # 1. Single Plate dCq
  output$download_plot_dCq <- create_qPCR_downloadHandler(
    input, 
    prepared_dCq_data, # Das neue Zwischen-Reactive!
    "dCq", 
    NULL, NULL, NULL   # Diese 3 sind jetzt NULL, da sie aus 'd' kommen
  )
  
  # 2. Single Plate ddCq
  output$download_plot_ddCq <- create_qPCR_downloadHandler(
    input = input, 
    data_reactive = prepared_ddCq_data, 
    type_suffix = "ddCq", 
    col_name = NULL, 
    err_name = NULL, 
    y_lab = NULL
  )
  
  # 3. Multi Plate (dCq)
  output$download_plot_MultiAbs <- create_qPCR_downloadHandler(
    input = input, 
    data_reactive = prepared_MultiAbs_data, 
    type_suffix = "Multi", # Prüfe ob deine Export-IDs im UI auch "plot_w_Multi" etc. heißen
    col_name = NULL, 
    err_name = NULL, 
    y_lab = NULL
  )
  
  # 4. Multi Plate Relative (ddCq)
  output$download_plot_MultiRel <- create_qPCR_downloadHandler(
    input = input, 
    data_reactive = prepared_MultiRel_data, 
    type_suffix = "DDCtMulti", # WICHTIG: Muss zu den UI-IDs passen (z.B. plot_w_DDCtMulti)
    col_name = NULL, 
    err_name = NULL, 
    y_lab = NULL
  )
  
  
  
}
