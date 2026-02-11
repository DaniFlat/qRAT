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
#' @importFrom plotly renderPlotly ggplotly plotlyOutput
#' @importFrom shinyWidgets virtualSelectInput
#' @importFrom shinyjqui orderInput


library("shiny")
library("bslib")
library("shinyWidgets")
library("HTqPCR")
library("ddCt")
library("plotly")
#library("shinyjs")
#library("ggplot2")
library("scales")
library("data.table")
library("DT")
library("dplyr")
library("waiter")
library("tidyr")
library("stringr")
library("magrittr")
library("shinycssloaders")
library("xfun")
library("viridisLite")
library("shinyjqui")
library("ggpubr")
library("rlang")

server <- function(input, output, session) {
  
 waiter_hide()
  
  
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
  
  
  # Marker, ob die aktuellen Dateien bereits validiert wurden
  MP_files_validated <- reactiveVal(FALSE)
  
  # Sobald neue Dateien ausgewählt werden, setzen wir den Marker zurück
  observeEvent(input$plates, {
    MP_files_validated(FALSE)
  })
  
  
  observeEvent(readMPData(), {
    req(readMPData())
    
    # Nur ausführen, wenn wir noch nicht validiert haben
    if (!MP_files_validated()) {
      ff <- input$plates
      filelist <- readMPData()$all_raw_list # Wir brauchen Zugriff auf die einzelnen DFs
      
      for (i in 1:nrow(ff)) {
        dt <- filelist[[i]]
        c_names <- tolower(names(dt))
        
        MP_validateInputFile(
          Sample_Column = length(intersect(sampleColumnList, c_names)) >= 1,
          Well_Column   = length(intersect(wellColumnList, c_names)) >= 1,
          Gene_Column   = length(intersect(geneColumnList, c_names)) >= 1,
          Cq_Column     = length(intersect(cqColumnList, c_names)) >= 1,
          plateName     = ff$name[i]
        )
      }
      # Jetzt auf TRUE setzen, damit bei Grouping-Änderungen Ruhe ist
      MP_files_validated(TRUE)
    }
  })
  
 
  
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
  
  
  
  
  
  
  # --- Hilfsfunktion für den Validierungs-Dialog (Multiple Plates) ---
  
  MP_validateInputFile <- function(Sample_Column, Well_Column, Gene_Column, Cq_Column, plateName) {
    
    # Sicherstellen, dass plateName ein einzelner String ist
    safePlateName <- if(length(plateName) > 0) as.character(plateName[1]) else "Unknown File"
    
    # Konvertierung zu logischen Einzelwerten
    is_sample <- isTRUE(as.logical(Sample_Column))
    is_well   <- isTRUE(as.logical(Well_Column))
    is_gene   <- isTRUE(as.logical(Gene_Column))
    is_cq     <- isTRUE(as.logical(Cq_Column))
    
    all_clear <- all(is_sample, is_well, is_gene, is_cq)
    
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
    
    showModal(modalDialog(
      title = tagList(icon("layer-group"), " Multi-Plate Validation"),
      size = "m",
      easyClose = FALSE,
      
      div(
        style = "padding: 10px;",
        div(
          style = "background: #325d88; color: white; padding: 12px 15px; border-radius: 8px 8px 0 0; display: flex; align-items: center; justify-content: space-between; border-bottom: 2px solid rgba(0,0,0,0.1);",
          tags$span(tags$b("Current File:"), style = "margin-right: 10px; opacity: 0.8;"),
          tags$code(safePlateName, style = "color: #ffca28; font-size: 1.1em; background: transparent; padding: 0;")
        ),
        
        if (all_clear) {
          div(class = "alert alert-success", style = "border-radius: 0 0 10px 10px; border: none; display: flex; align-items: center; margin-bottom: 20px;",
              icon("circle-check", style = "margin-right: 15px; font-size: 1.5rem;"),
              div(tags$b("Validation Successful"), br(), "All required columns were found.")
          )
        } else {
          div(class = "alert alert-danger", style = "border-radius: 0 0 10px 10px; border: none; display: flex; align-items: center; margin-bottom: 20px;",
              icon("triangle-exclamation", style = "margin-right: 15px; font-size: 1.5rem;"),
              div(tags$b("Columns Missing!"), br(), "Please check the column headers or delimiter.")
          )
        },
        
        div(
          style = "background: #f8f9fa; border-radius: 10px; border: 1px solid #e9ecef; overflow: hidden;",
          status_row("Sample Column", is_sample, "vial"),
          status_row("Well Column", is_well, "border-all"),
          status_row("Gene / Target Column", is_gene, "dna"),
          status_row("Cq / Ct Column", is_cq, "chart-line")
        )
      ),
      
      footer = tagList(
        modalButton("Next / Close", icon = icon("chevron-right"))
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
  
  # Single Plate read file
  
  readSinglePlateData <- reactive({
    
    if (input$Sep == "tab") sep <- "\t"
    if (input$Sep == "comma") sep <- ","
    if (input$Sep == "semicolon") sep <- ";"
    if (input$Sep == "auto") sep <- "auto"
    
    req(input$dtfile)
    f <- input$dtfile
    if (is.null(f)) {
      return(NULL)
    } else {
      filex <- f$datapath
    }
    
    #Check file type
    req(input$dtfile)
    ext <- file_ext(f$datapath)
    if (!ext %in% c(
      'text/csv',
      'text/comma-separated-values',
      'text/tab-separated-values',
      'text/plain',
      'csv',
      'txt',
      'tsv')) {
      showModal(check_filetype(failed = TRUE))
    }
    
    # Scan through beginning of file, max 100 lines and look for row with "Well"
    file.header <- readLines(con = filex, n = 100)
    n.header <- grep("^Well", file.header) - 1
    if (length(n.header) == 0) {
      n.header <- 0
    }
    
    SP_data <- fread(filex, header = TRUE, sep = sep, skip = n.header)
    
    
    ##Validate Input Data File
    #change all column names to lowercase for easy comparison with prepared column-name-lists
    SP_data <- SP_data %>%
      rename_with(tolower)
    
    column_names <- colnames(SP_data)
    
    validateInputFile <- function(Sample_Column=TRUE, Well_Column=TRUE, Gene_Column=TRUE, Cq_Column=TRUE) {
      if (length(intersect(sampleColumnList,column_names)) < 1) Sample_Column = FALSE
      if (length(intersect(wellColumnList,column_names)) < 1) Well_Column = FALSE
      if (length(intersect(geneColumnList,column_names)) < 1) Gene_Column = FALSE
      if (length(intersect(cqColumnList,column_names)) < 1) Cq_Column = FALSE
      
      showModal(SP_validateInputFile(Sample_Column, Well_Column, Gene_Column, Cq_Column))
    }
    
    validateInputFile()
    
    
    ## change all column names to lowercase and rename columns to Well, Sample, Gene, Ct
    SP_data <- SP_data %>%
      rename(
        Sample = any_of(sampleColumnList),
        Well = any_of(wellColumnList),
        Gene = any_of(geneColumnList),
        Ct = any_of(cqColumnList)
      )
    
    
    ## generate replicate number for dataview
    SP_data$TempRepNum <- paste(SP_data$Sample, SP_data$Gene)
    SP_data$rp.num <- ave(SP_data$Sample, SP_data$TempRepNum, FUN = seq_along)
    SP_data <- subset(SP_data, select = c(Well, Sample, Gene, Ct, rp.num)) %>%
      filter(Well != "")
    
    ## replace "," with "." for dataview
    CorrectCt <- SP_data$Ct
    CorrectCt <- gsub(",", ".", CorrectCt, ignore.case = TRUE)
    SP_data$Ct <- as.numeric(CorrectCt)
    
    # Update Inputfield for NTC Selection with all Samplenames
    NTCs <- unique(as.character(SP_data$Sample))
    updateVirtualSelect("NTC_Input", choices = NTCs)
    
    # raw data for plots (spatial plate view; cq distribution)
    plateView <- SP_data %>%
      separate(Well,
               into = c("text", "num"),
               sep = "(?<=[A-Za-z])(?=[0-9])"
      )
    
    SP_data
    
    return(list(data = SP_data, PV = plateView))
  })
  
  
  
  setData <- reactive({
    req(input$dtfile)
    info <- readSinglePlateData()
    xx <- info$data
    if (is.null(xx)) return(NULL)
    
    sep_val <- if (input$Sep == "tab") "\t" else if (input$Sep == "comma") "," else if (input$Sep == "semicolon") ";" else "\t"
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
      processed <- processed %>%
        mutate(Ct = remove_replicates(Ct)) %>%
        mutate(Ct = remove_maxCt(Ct))
      
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
    if (isTRUE(input$groupSamples) && !is.null(confirmedGrouping())) {
      groups <- confirmedGrouping()
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
    
    # Daten schreiben für das qPCR-Paket
    fwrite(processed, file = "fileSingle.csv", sep = sep_val)
    
    # Downstream Updates
    Genes <- unique(processed$Gene)
    refs <- grep("^(ACT|UBQ|REF|RF|BACTIN|ACTIN|GAPDH|TUB|TUBULIN).*", Genes, ignore.case = TRUE)
    selected_refs <- if (length(refs) == 0) Genes[length(Genes)] else Genes[refs]
    updateVirtualSelect("Refs", choices = Genes, selected = selected_refs)
    
    # qPCR Objekt Erstellung
    htset <- NULL; ddset <- NULL; Samples <- NULL
    xdata <- read.qPCRtable("fileSingle.csv", sep = sep_val)
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))
      updateVirtualSelect("Mock", choices = Samples, selected = Samples[1])
      updateVirtualSelect("Comps", choices = Samples, selected = Samples)
    }
    
    return(list(
      data = processed, 
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
  remove_replicates <- function(x, na.rm = TRUE, ...) {
    repDev <- input$repDev
    
    val <- mean(x, na.rm = na.rm)
    y <- x
    y[x < (val - repDev)] <- NA
    y[x > (val + repDev)] <- NA
    y
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
    repDev <- input$repDevMulti
    
    val <- mean(x, na.rm = na.rm)
    y <- x
    y[x < (val - repDev)] <- NA
    y[x > (val + repDev)] <- NA
    y
  }
  
  
  
  # --- Einlese-Funktion (Multi Plate) ---
  readMPData <- reactive({
    req(input$plates)
    
    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"
    if (input$SepM == "auto") sep <- "auto"
    
    ff <- input$plates
    numfiles <- nrow(ff)
    filelist <- list()
    
    # 1. Dateityp-Check
    for (i in 1:numfiles) {
      ext <- tools::file_ext(ff$name[i])
      if (!ext %in% c("csv", "txt", "tsv")) {
        showModal(modalDialog(
          title = "Invalid File Type",
          paste("File", ff$name[i], "is not a supported format (.csv, .txt, .tsv)"),
          easyClose = TRUE
        ))
        return(NULL)
      }
    }
    
    # 2. Dateien einzeln einlesen und validieren
    for (i in 1:numfiles) {
      file_path <- ff$datapath[i]
      current_filename <- ff$name[i]
      
      # Header-Erkennung
      file.header <- readLines(con = file_path, n = 100)
      n.header <- grep("^Well", file.header) - 1
      if (length(n.header) == 0) n.header <- 0
      
      # Einlesen
      dt <- fread(file_path, header = TRUE, sep = sep, skip = n.header)
      filelist[[i]] <- dt
      
      # # Validierung Trigger
      # c_names <- tolower(names(dt))
      # MP_validateInputFile(
      #   Sample_Column = length(intersect(sampleColumnList, c_names)) >= 1,
      #   Well_Column   = length(intersect(wellColumnList, c_names)) >= 1,
      #   Gene_Column   = length(intersect(geneColumnList, c_names)) >= 1,
      #   Cq_Column     = length(intersect(cqColumnList, c_names)) >= 1,
      #   plateName     = current_filename
      # )
    }
    
    # 3. Datenverarbeitung
    multiplePlates <- rbindlist(filelist, use.names = TRUE, fill = TRUE, idcol = "PlateNumber")
    
    multiplePlates <- multiplePlates %>%
      rename_with(tolower) %>%
      rename(
        Sample = any_of(sampleColumnList),
        Well = any_of(wellColumnList),
        Gene = any_of(geneColumnList),
        Ct = any_of(cqColumnList),
        PlateNumber = platenumber
      )
    
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
    MultiplePlatesView <- multiplePlates %>%
      separate(Well,
               into = c("text", "num"),
               sep = "(?<=[A-Za-z])(?=[0-9])"
      )
    
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
        mutate(Ct = remove_replicatesMP(Ct)) %>% 
        mutate(Ct = remove_maxCtMulti(Ct))
      
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
    updateVirtualSelect("IPC", choices = SamplesWithIPCs)
    
    # Spalte OriginalRowID entfernen, falls sie später stört
    correctRep$OriginalRowID <- NULL
    
    return(list(data = correctRep, badReplicates = badRep))
  })
  
  # --- multiData 2: Trennung der IPCs von den Haupt-Proben ---
  multiData2 <- reactive({
    req(input$plates)
    
    sep <- if (input$SepM == "comma") "," else if (input$SepM == "semicolon") ";" else "\t"
    
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
    
    fwrite(info2, file = "fileMultiNoIPCs.csv", sep = sep)
    
    # Gene einlesen und Housekeeper vorschlagen
    Genes <- unique(info2$Gene)
    refs <- grep("^(ACT|UBQ|REF|RF|BACTIN|ACTIN|GAPDH|TUB|TUBULIN).*", Genes, ignore.case = TRUE)
    selected_refs <- if (length(refs) == 0) Genes[length(Genes)] else Genes[refs]
    updateVirtualSelect("RefsM", choices = Genes, selected = selected_refs)
    
    # Konvertierung in das interne qPCR-Format
    htset <- NULL; ddset <- NULL; Samples <- NULL
    xdata <- read.qPCRtableMulti("fileMultiNoIPCs.csv", sep = sep)
    
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))
      
      updateVirtualSelect("MockM", choices = Samples, selected = Samples[1])
      updateVirtualSelect("CompsM", choices = Samples, selected = Samples)
      
      # Dynamische Anpassung der Plot-Achsen
      if (length(Samples) > 4 | max(nchar(Samples)) > 6) {
        updateSliderInput(session, "side1", value = 8)
        updateSliderInput(session, "xsrt", value = 30)
      }
    }
    
    return(list(
      SamplesNoIPCs = info2, 
      xdata = xdata, 
      data.ht = htset, 
      data.ddct = ddset, 
      CompSamples = Samples
    ))
  })
  
  
  
  
  ## HTqPCR Single Plate
  resultHTqPCR <- reactive({
    info <- setData()
    
    if (is.null(info)) {
      return(NULL)
    }
    
    if (input$limma_input == "dCq") {
      dataset <- info$data.ht
    } else {
    }
    
    info2 <- info$data.ht
    Samples <- as.character(info2$Sample)
    refs <- input$Refs
    multiComps <- listOfCompsMultiple()
    adjustMethod <- input$adjustMethod
    
    if (input$compType == "Single Comparison") {
      comp.type <- 1
      comps <- input$Comps
    } else {
      comp.type <- 2
      comps <- multiComps
    }
    
    if (input$adjustMethod == "Benjamini & Hochberg") {
      adjustMethod <- "BH"
    } else {
      if (input$adjustMethod == "Holm") {
        adjustMethod <- "holm"
      }
      adjustMethod <- "bonferroni"
    }
    
    res <- calHTqPCR(dataset, refs, comps, comp.type, adjustMethod)
    return(res)
  })
  
  ## HTqPCR Multiple Plates
  resultHTqPCRMulti <- reactive({
    info <- multiData2()
    if (is.null(info)) {
      return(NULL)
    }
    
    if (input$limma_inputMulti == "dCq") {
      dataset <- info$data.ht
    } else {
    }
    
    info2 <- info$data.ht
    Samples <- as.character(info2$Sample)
    refs <- input$RefsM
    multiComps <- listOfCompsMultipleM()
    
    if (input$compTypeM == "Single Comparison") {
      comp.type <- 1
      comps <- input$CompsM
    } else {
      comp.type <- 2
      comps <- multiComps
    }
    
    if (input$adjustMethodMulti == "Benjamini & Hochberg") {
      adjustMethod <- "BH"
    } else {
      if (input$adjustMethod == "Holm") {
        adjustMethod <- "holm"
      }
      adjustMethod <- "bonferroni"
    }
    
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
      filter(FC != "") %>%
      filter(FC.sd != "")
    
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
      filter(expr != "") %>%
      filter(expr.sd != "")
    expr.abs$"-dCq" <- expr.abs$"-dCq" * -1
    SamplesdCt <- unique(as.character(expr.abs$Sample))
    SamplesddCt <- unique(as.character(expr.rel$Sample))
    GenesdCt <- unique(as.character(expr.abs$Gene))
    GenesddCt <- unique(as.character(expr.rel$Gene))
    
    # Update Inputfields (Plot Settings and Input Settings)
    updateVirtualSelect("SamplePicker", choices = SamplesdCt, selected = SamplesdCt)
    updateOrderInput(session, "x_order", items = SamplesdCt)
    updateVirtualSelect("SamplePickerDDCt", choices = SamplesddCt, selected = SamplesddCt)
    updateOrderInput(session, "x_orderDDCt", items = SamplesddCt)
    updateVirtualSelect("GenePicker", choices = GenesdCt, selected = GenesdCt)
    updateVirtualSelect("GenePickerDDCt", choices = GenesddCt, selected = GenesddCt)
    
    #prepare ddCt output for statistical analysis
    ddCt_stat <- expr.rel %>% select(Sample, Gene, ddCt)
    ddCt_stat <- rename(ddCt_stat, Ct = ddCt)
    
    
    fwrite(ddCt_stat, file = "ddCtsingle.csv")
    
    ## read to qPCR format
    ddCt_stat_set <- NULL
    xdata <- read.qPCRtable("ddCtsingle.csv")
    if (!is.null(xdata)) {
      ddCt_stat_set <- xdata$data.ht
    }
    
    
    list(relative = expr.rel, absolute = expr.abs, stat_set = ddCt_stat_set)
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
      filter(FC != "") %>%
      filter(FC.sd != "")
    
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
      filter(expr != "") %>%
      filter(expr.sd != "")
    
    expr.abs$"-dCq" <- expr.abs$"-dCq" * -1
    SamplesdCtMulti <- unique(as.character(expr.abs$Sample))
    SamplesddCtMulti <- unique(as.character(expr.rel$Sample))
    
    GenesdCtMulti <- unique(as.character(expr.abs$Gene))
    GenesddCtMulti <- unique(as.character(expr.rel$Gene))
    
    # Update Inputfields (Plot Settings and Input Settings)
    updateVirtualSelect("SamplePickerMulti", choices = SamplesdCtMulti, selected = SamplesdCtMulti)
    updateOrderInput(session, "x_orderMulti", items = SamplesdCtMulti)
    updateVirtualSelect("SamplePickerDDCtMulti", choices = SamplesddCtMulti, selected = SamplesddCtMulti)
    updateOrderInput(session, "x_orderDDCtMulti", items = SamplesddCtMulti)
    updateVirtualSelect("GenePickerMulti", choices = GenesdCtMulti, selected = GenesdCtMulti)
    updateVirtualSelect("GenePickerDDCtMulti", choices = GenesddCtMulti, selected = GenesddCtMulti)
    
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
    req(multiData2())
    x <- multiData2()$SamplesNoIPCs
    
    if (isTRUE(input$Id027)) {
      calF <- extractIPCs()
      req(calF$CalibrationFactors)
      
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
    
    data_orig <- multiData2()
    req(data_orig)
    noncalibrated_df <- data_orig$SamplesNoIPCs
    
    # Hilfsspalte für die Auswahl erstellen
    cal <- calibrated_df %>% unite("Sample_Gene", Sample:Gene, remove = FALSE)
    non <- noncalibrated_df %>% unite("Sample_Gene", Sample:Gene, remove = FALSE)
    
    return(list(calibrated = cal, noncalibrated = non))
  })
  
  # Separater Observer für das Update des Dropdowns
  observe({
    info <- interPlateCalibration_Plot()
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
    if (is.null(info)) NULL
    
    info$data  %>% rename(Replicate = rp.num)
    
  })
  
  ## Raw Data Table Plate Plot Spatial Original Data
  
  output$plotCtCard <- renderPlotly({
    info <- readSinglePlateData()
    hovertext <- info$PV
    
    if (is.null(info)) {
      NULL
    } else {
      fig <- plot_ly(as.data.frame(info$PV), x = ~num, y = ~text, z = ~Ct, type = "heatmap", xgap = 1, ygap = 1, colors = "RdYlBu", colorbar = list(title = "Cq"),
                     hoverinfo = 'text',
                     text = ~paste("Sample:", hovertext$Sample,
                                   "<br>Gene:", hovertext$Gene,
                                   "<br>Cq:", hovertext$Ct)) %>%
        layout(
          title = list(text = "Plate View", y=1),
          font=t,
          xaxis = list(side = "top", title = "", type = "linear", tickmode = "linear", showgrid = FALSE, ticks = "outside", tickwidth=2),
          yaxis = list(autorange = "reversed", title = "", showgrid = FALSE, ticks = "outside", tickwidth=2)
        ) %>%
        config(
          displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
          toImageButtonOptions = list(format = "svg", scale = 3)
        )
    }
  })
  
  ## Raw Data Table Cq distribution Original Data
  
  output$plotCtDistrib <- renderPlotly({
    info <- readSinglePlateData()
    
    if (is.null(info)) NULL
    
    fig <- ggplot(info$PV, aes(x = Ct, colour = Sample), size=0.6, scale_color_viridis()) +
      geom_density() +
      geom_rug() +
      scale_color_brewer(palette = "Set1")+
      #theme_pubr(base_size=4.76 * .pt)+
      list(ggplottheme)
    
    ggplotly(fig) %>%
      layout(
        title = list(text = "Cq Distribution", y=1),
        plot_bgcolor = "transparent",
        font=t,
        yaxis = list(title = "Density", linewidth = 2, ticks = "outside", tickwidth=2),
        xaxis = list(title = "Cq", linewidth = 2, ticks = "outside", tickwidth=2),
        legend=list(title=list(text='Samples'))
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "png", scale = 3), scrollZoom = TRUE
      )
  })
  
  
  ## Single Plate Show bad replicates Table
  
  output$dataSinglePlateBadRep <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- setData()
    data <- info$badReplicates %>%
      select(Well, Sample, Gene, rp.num, Cq) %>%
      rename(Replicate = rp.num)
    
    if (is.null(info)) NULL else data
  })
  
  
  ## Single Plate Filter Plot (how many replicates are neglected per sample)
  
  output$SinglePlateFilterPlot <- renderPlotly({
    info <- setData()
    if (is.null(info)) NULL
    
    filtered <- info$badReplicates %>%
      count(Sample)
    
    unfiltered <- info$data %>%
      filter(!is.na(Ct)) %>%
      count(Sample)
    
    together <- bind_rows(filtered, unfiltered, .id = "status") %>%
      mutate(status = replace(status, status == 1, "unreliable")) %>%
      mutate(status = replace(status, status == 2, "OK")) %>%
      plot_ly(
        x = ~Sample, y = ~n, type = "bar",
        name = ~status, color = ~status, colors = "Paired", marker = list(line = list(color = "rgba(0, 0, 0, .8)", width = 2))
      ) %>%
      layout(
        title = list(text = "filtered/unfiltered Replicates [n]", y=1),
        font=t,
        yaxis = list(title = "Count", showgrid = FALSE, showline = TRUE, zeroline = FALSE, ticks = "outside", linewidth = 2),
        barmode = "stack",
        xaxis = list(linewidth = 2, zeroline = FALSE, ticks = "outside", tickvals = ~Sample, ticktext = ~Sample)
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "svg", scale = 3)
      )
  })
  
  
  
  ## Single Plate Boxplot after filtering
  
  output$SinglePlateBoxplot <- renderPlotly({
    info <- setData()
    if (is.null(info)) NULL
    
    finalData <- as.data.frame(info$data) %>%
      plot_ly(
        x = ~Sample, y = ~Ct, color = ~Gene, type = "box", colors = "Spectral", boxpoints = "all", jitter = 0.3, pointpos = 0,
        marker = list(line = list(color = "rgba(0, 0, 0, .8)", width = 1)),
        line = list(color = "rgba(0, 0, 0, .8)", width = 1)
      ) %>%
      layout(
        title = list(text = "Variation within Samples", y=1),
        legend = list(title = list(text = "Genes")),
        font=t,
        yaxis = list(title = "Cq", showgrid = FALSE, showline = TRUE, linewidth = 2, mirror = TRUE, ticks = "outside"),
        boxmode = "group",
        xaxis = list(showgrid = FALSE, showline = TRUE, linewidth = 2, mirror = TRUE, ticks = "outside")
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "svg", scale = 3), scrollZoom = TRUE, displayModeBar = TRUE
      )
  })
  
  ## dCt Expression Data Table Single Plate
  
  output$ddctAbsolute <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdct()
    if (is.null(info)) {
      return(NULL)
    }
    
    info$absolute %>%
      mutate(across(where(is.numeric), round, 2)) %>%
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
    
    formatChoice   <- input$exportFormat %||% "png"
    PlotWidth      <- get_input("width") %||% 700
    PlotHeight     <- get_input("height") %||% 500
    PlotScale      <- get_input("scale_factor") %||% 1
    
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
      colourpalette <- viridis::viridis(numberOfGenes)
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
      "<br><b>Error:</b> ", round(err_vec, 3)
    )
    
    # 6. Plot-Basis erstellen
    p <- ggplot(df, aes(x = Sample, y = !!sym(PlotDataCol), fill = Gene, color = Gene, group = Gene, text = tooltip_text))
    
    if (PlotType == "Bar Chart") {
      p <- p + geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7)
    } else if (PlotType == "Boxplot") {
      p <- p + geom_boxplot(position = position_dodge(0.9))
    } else {
      p <- p + geom_point(position = position_dodge(0.9), size = 3)
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
        show.legend = FALSE
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
      stats_df$signif <- symnum(stats_df$adj.p.value, corr = FALSE, na = FALSE, 
                                cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                symbols = c("***", "**", "*", "ns"))
      
      df_merged <- merge(df, stats_df[, c("genes", "Label", "signif")], 
                         by.x = c("Gene", "Sample"), by.y = c("genes", "Label"), all.x = TRUE)
      
      p <- p + geom_text(
        data = df_merged,
        aes(x = Sample, y = !!sym(PlotDataCol) + !!sym(PlotDataError), label = signif, group = Gene),
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
    
    # 11. Plotly Export
    if (isTRUE(export_mode)) {
      return(p) # Gibt das ggplot-Objekt zurück
    } else {
      return(
        ggplotly(p, tooltip = "text") %>% 
          layout(margin = list(l = 70, r = 50, b = 100, t = 60)) %>%
          config(displaylogo = FALSE)
      )
    }
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
    return(list(
      plot_df = df_plot,
      stats = calLimma(),
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  })
  
  
  output$dcqPlot <- renderPlotly({
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
    
    # Rückgabe als Liste für Plot und Export
    return(list(
      plot_df = df_plot,
      stats = calLimma(),
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  }) 
  
  
 
  
  output$ddctRelGraph2 <- renderPlotly({
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
      mutate(across(where(is.numeric), round, 2)) %>%
      rename(ddCq = ddCt, ddCq.sd = ddCt.sd)
  })
  
  
  
  
  
  ## MULTIPLE PLATES OUTPUT
  
  
  ## Raw Data Table Multiple Plates
  
  output$multiplePlatesData <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- readMPData()
    if (is.null(info)) NULL else info$data %>% rename(Plate = PlateNumber, Replicate = rp.num)
  })
  
  
  ## Raw Data Table Plate Plot Spatial Original Data Multiple Plates
  
  output$MultiplotCtCard <- renderPlotly({
    info <- readMPData()
    if (is.null(info)) {
      NULL
    } else {
      plateInput <- input$PlateSelect
    }
    
    plateData <- info$MPV %>% filter(PlateNumber == plateInput)
    
    plot_ly(plateData, x = ~num, y = ~text, z = ~Ct, type = "heatmap", colors = "RdYlBu", xgap = 1, ygap = 1, colorbar = list(title = "Cq"),
            hoverinfo = 'text',
            text = ~paste("Sample:", plateData$Sample,
                          "<br>Gene:", plateData$Gene,
                          "<br>Cq:", plateData$Ct)) %>%
      layout(
        title = list(text = "Plate View", y=1),
        font=t,
        xaxis = list(side = "top", title = "", type = "linear", tickmode = "linear", showgrid = FALSE, ticks = "outside", tickwidth=2),
        yaxis = list(autorange = "reversed", title = "", showgrid = FALSE, ticks = "outside", tickwidth=2)
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "svg", scale = 3)
      )
  })
  
  ## Raw Data Table Cq distribution Original Data Multiple Plates
  
  output$MultiplotCtDistrib <- renderPlotly({
    info <- multiData()
    
    if (is.null(info)) NULL
    
    fig <- ggplot(info$data, aes(x = Ct, colour = Sample), size=0.6, scale_color_viridis()) +
      geom_density() +
      geom_rug() +
      scale_color_brewer(palette = "Set1")+
      #theme_pubr(base_size=4.76 * .pt)+
      list(ggplottheme)
    
    ggplotly(fig) %>%
      layout(
        title = list(text = "Cq Distribution", y=1),
        plot_bgcolor = "transparent",
        font=t,
        yaxis = list(title = "Density", linewidth = 2, ticks = "outside", tickwidth=2),
        xaxis = list(title = "Cq", linewidth = 2, ticks = "outside", tickwidth=2),
        legend=list(title=list(text='Samples'))
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "svg", scale = 3), scrollZoom = TRUE
      )
  })
  
  
  
  ## Data Table Multiple Plates Show bad replicates
  
  output$dataMultiplePlatesBadRep <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- multiData()
    data <- info$badReplicates %>%
      select(PlateNumber, Well, Sample, Gene, rp.num, Cq) %>%
      rename(Plate = PlateNumber, Replicate = rp.num)
    if (is.null(info)) NULL else data
  })
  
  
  
  ## Multiple Plates Filter Plot (how many replicates are neglected per sample)
  
  output$MultiplePlatesFilterPlot <- renderPlotly({
    info <- multiData()
    if (is.null(info)) NULL
    
    filtered <- info$badReplicates %>%
      count(Sample)
    
    unfiltered <- info$data %>%
      filter(!is.na(Ct)) %>%
      count(Sample)
    
    together <- bind_rows(filtered, unfiltered, .id = "status") %>%
      mutate(status = replace(status, status == 1, "unreliable")) %>%
      mutate(status = replace(status, status == 2, "OK")) %>%
      plot_ly(
        x = ~Sample, y = ~n, type = "bar",
        name = ~status, color = ~status, colors = "Paired", marker = list(line = list(color = "rgba(0, 0, 0, .8)", width = 2))
      ) %>%
      layout(
        title = list(text = "filtered/unfiltered Replicates [n]", y=1),
        font=t,
        yaxis = list(title = "Count", showgrid = FALSE, showline = TRUE, zeroline = FALSE, ticks = "outside"),
        barmode = "stack",
        xaxis = list(linewidth = 2, zeroline = FALSE, ticks = "outside", tickvals = ~Sample, ticktext = ~Sample)
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "svg", scale = 3)
      )
  })
  
  ## Multiple Plate Boxplot after filtering
  
  output$MultiplePlatesBoxplot <- renderPlotly({
    info <- multiData()
    if (is.null(info)) NULL
    
    finalData <- as.data.frame(info$data) %>%
      plot_ly(
        x = ~Sample, y = ~Ct, color = ~Gene, type = "box", colors = "Spectral", boxpoints = "all", jitter = 0.3, pointpos = 0,
        marker = list(line = list(color = "rgba(0, 0, 0, .8)", width = 1)),
        line = list(color = "rgba(0, 0, 0, .8)", width = 1)
      ) %>%
      layout(
        title = list(text = "Variation within Samples", y=1),
        legend = list(title = list(text = "Genes")),
        font=t,
        yaxis = list(title = "Cq", showgrid = FALSE, showline = TRUE, linewidth = 2, mirror = TRUE, ticks = "outside"),
        boxmode = "group",
        xaxis = list(showgrid = FALSE, showline = TRUE, linewidth = 2, mirror = TRUE, ticks = "outside")
      ) %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = "svg", scale = 3), scrollZoom = TRUE, displayModeBar = TRUE
      )
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
      data <- info$CalibrationFactors %>% mutate(across(where(is.numeric), round, 4))
      if (is.null(info)) NULL else data
    } else {}
  })
  
  ## Inter-Plate Ct values, calibrated with calibration factor
  
  output$interPlateCalibration <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    if (input$Id027 == TRUE) {
      info <- interPlateCalibration()
      data <- info %>%
        select(PlateNumber, Well, Sample, Gene, Ct, rp.num) %>%
        mutate(across(where(is.numeric), round, 2)) %>%
        rename(Plate = PlateNumber, Cq = Ct, Replicate = rp.num)
      if (is.null(info)) NULL else data
    } else {}
  })
  
  
 
  
  
  
  #IPC Comparison Plot
  output$plotIPCcomparison <- renderPlotly({
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
      geom_boxplot(aes(fill = Status), alpha = 0.1, outlier.shape = NA, color = "grey70") +
      # Einzelne Punkte (Replikate)
      geom_jitter(width = 0.15, size = 2, alpha = 0.8, 
                  aes(text = paste("Plate:", PlateNumber, "<br>Ct:", round(Ct, 2), "<br>Gene:", Gene))) +
      
      # scales = "free_y" erlaubt jedem Gen seine eigene Y-Achse
      facet_wrap(~Sample_Gene, scales = "free_y") +
      theme_minimal() +
      scale_color_brewer(palette = "Set1") + # Bessere Farben für Platten
      labs(title = "Effect of Inter-Plate Calibration",
           x = "", y = "Ct Value", color = "Plate") +
      theme(strip.background = element_rect(fill = "#f8f9fa"),
            legend.position = "bottom")
    
    # Umwandlung in interaktiven Plot
    ggplotly(p, tooltip = "text") %>% 
      layout(margin = list(l = 50, r = 50, b = 50, t = 80))
  })
  
  
  
  
  
  ## dCt Gene Expression Table Multiple Plates
  
  output$ddctAbsoluteMulti <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdctMulti()
    if (is.null(info)) {
      return(NULL)
    }
    info$multiAbsolute %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      rename(dCq = dCt, dCq.sd = dCt.sd, RQ = expr, RQ.sd = expr.sd)
  })
  
 
  
  ## ddCq Gene Expression Table Multiple Plates
  
  output$ddctRelativeMulti <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- resultDdctMulti()
    if (is.null(info)) {
      return(NULL)
    }
    info$multiRelative %>%
      mutate(across(where(is.numeric), round, 2)) %>%
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
    
    # Rückgabe
    return(list(
      plot_df = df_plot,
      stats = calLimmaMulti(),
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  })
  
  
  output$ddctMultiAbsGraph <- renderPlotly({
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
    
    return(list(
      plot_df = df_plot,
      stats = calLimmaMulti(),
      col_val = col_val,
      err_col = err_col,
      ylab = ylab
    ))
  })  
  

  output$ddctRelGraphMulti <- renderPlotly({
    d <- prepared_MultiRel_data()
    # Nutzt input_type "MultiRel" -> Suffix "DDCtMulti"
    render_qPCR_plot(d$plot_df, d$stats, d$col_val, d$err_col, d$ylab, input_type = "MultiRel")
  })  
  
  
  
  ## MIQE Checks
  ## Single Plate MIQE check
  
  output$SinglePlateMIQEcheck <- renderUI({
    info <- setData()
    
    if (is.null(info)) {
      NULL
    } else {
      countReplicates <- info$data %>%
        filter(!is.na(Ct)) %>%
        count(Sample, Gene) %>%
        filter(n < 2)
      
      print(countReplicates)
    }
    
    output$countTable <- renderDT(countReplicates)
    
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
  
  output$MultiplePlateMIQEcheck <- renderUI({
    info <- multiData()
    
    if (is.null(info)) {
      NULL
    } else {
      countReplicates <- info$data %>%
        filter(!is.na(Ct)) %>%
        count(Sample, Gene) %>%
        filter(n < 2)
    }
    
    output$countTable <- renderDT(countReplicates)
    
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
  
  
  
  ## Limma Outputs
  
  output$resultLimma <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    dt <- calLimma()
    if (is.null(dt)) {
      return(NULL)
    }
    
    dt %>%
      select(Label, genes, t.test, p.value, adj.p.value) %>%
      rename(Sample = Label, Gene = genes) %>%
      mutate(significance = case_when(
        adj.p.value > 0.05 ~ "",
        adj.p.value > 0.01 ~ "*",
        adj.p.value > 0.001 ~ "**",
        adj.p.value > 0.0001 ~ "***",
        adj.p.value > 0 ~ "****")) %>%
      mutate(across(c("p.value", "adj.p.value"), label_pvalue(accuracy = 0.0001))) %>%
      mutate(across(where(is.numeric), round, 4))
  })
  
  output$resultLimmaMulti <- renderDT(extensions = c("Buttons"), options = table_options(), rownames = FALSE, filter = "top", {
    dt <- calLimmaMulti()
    if (is.null(dt)) {
      return(NULL)
    }
    dt %>%
      select(Label, genes, t.test, p.value, adj.p.value) %>%
      rename(Sample = Label, Gene = genes) %>%
      mutate(significance = case_when(
        adj.p.value > 0.05 ~ "",
        adj.p.value > 0.01 ~ "*",
        adj.p.value > 0.001 ~ "**",
        adj.p.value > 0.0001 ~ "***",
        adj.p.value > 0 ~ "****")) %>%
      mutate(across(c("p.value", "adj.p.value"), label_pvalue(accuracy = 0.0001))) %>%
      mutate(across(where(is.numeric), round, 4))
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
  
  listOfCompsMultiple <- reactive({
    n <- counter$n
    comparisonsList <- NULL
    if (n > 0) {
      isolate({
        multipleComparisons <- sapply(seq_len(n), function(i) {
          comparisonsList <- toString(AllInputs()[[paste0("CompsMultiple", i)]])
          str_replace(comparisonsList, ", ", ":")
        })
        multiComps <- paste(multipleComparisons, collapse = ",")
        multiComps <- gsub(" ", ".", multiComps, ignore.case = TRUE)
      })
    } else {
      NULL
    }
  })
  
  output$comparisonBoxes_ui <- renderUI({
    comparisonBoxes()
  })
  
  ### Comparison Boxes for Statistical Analysis
  ### Multiple Plates
  # Track the number of input boxes to render
  counterM <- reactiveValues(n = 0)
  
  # Track all user inputs
  AllInputs <- reactive({
    x <- reactiveValuesToList(input)
  })
  
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
    n <- counterM$n
    comparisonsListM <- NULL
    if (n > 0) {
      isolate({
        multipleComparisonsM <- sapply(seq_len(n), function(i) {
          comparisonsListM <- toString(AllInputs()[[paste0("CompsMultipleM", i)]])
          str_replace(comparisonsListM, ", ", ":")
        })
        multiCompsM <- paste(multipleComparisonsM, collapse = ",")
        multiCompsM <- gsub(" ", ".", multiCompsM, ignore.case = TRUE)
      })
    } else {
      NULL
    }
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
    info <- readSinglePlateData()
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
    info <- readMPData()
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
  
  observe({
    for (id in active_group_ids()) {
      observeEvent(input[[paste0("remove_grp_btn_", id)]], {
        removeUI(selector = paste0("#group_container_", id))
        active_group_ids(active_group_ids()[active_group_ids() != id])
      }, ignoreInit = TRUE, once = TRUE)
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
  observeEvent(input$open_grouping_modal_multi, {
    
    req(multiData()) 
  
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
    if (!is.null(multiData())) {
      # Wir greifen auf die Spalte 'Sample' im Daten-Element der Liste zu
      available_samples <- sort(unique(as.character(multiData()$data$Sample)))
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
    if (!is.null(multiData())) {
      available_samples <- sort(unique(as.character(multiData()$data$Sample)))
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
  
  # 4. Löschen von Gruppen (dynamische Observer)
  observe({
    for (id in active_group_ids_multi()) {
      observeEvent(input[[paste0("remove_grp_btn_multi_", id)]], {
        removeUI(selector = paste0("#group_container_multi_", id))
        active_group_ids_multi(active_group_ids_multi()[active_group_ids_multi() != id])
      }, ignoreInit = TRUE, once = TRUE)
    }
  })
  
  

  #Reference Gene Finder
  # 1. Dynamische Auswahl der Reference-Genes Single Plate
  observe({
    # Wir nehmen die Daten nach dem Filtering
    req(setData())
    available_genes <- sort(unique(as.character(setData()$data$Gene)))
    
    updateVirtualSelect("ref_candidates", choices = available_genes)
  })
  
  # Dynamische Auswahl der Reference-Genes Multi Plate
  observe({
    req(multiData2())
    # Wir nehmen die Gene aus dem bereinigten Multi-Plate Datensatz
    available_genes <- sort(unique(as.character(multiData2()$SamplesNoIPCs$Gene)))
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
      res <- run_normfinder_logic(ct_mat)
      setProgress(1)
    })
    
    return(res)
  })
  
  # Plot für die Stabilität
  output$ref_stability_plot <- renderPlotly({
    res <- validationResults()
    req(res)
    
    p <- ggplot(res, aes(x = reorder(Gene, Stability), y = Stability, fill = Stability)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      scale_fill_gradient(low = "#2ecc71", high = "#e74c3c") +
      theme_minimal() +
      labs(x = "Candidate Genes", y = "Stability Value (lower is better)",
           title = "NormFinder Ranking") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
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
      formatRound("Stability", 3)
  }) 
  
  
  ##Reference Gene Finder Multiplate
  validationResultsMulti <- eventReactive(input$run_val_btn_multi, {
    req(multiData2())
    candidates <- input$ref_candidates_multi
    req(candidates)
    
    raw_df <- multiData2()$SamplesNoIPCs
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
      res <- run_normfinder_logic(ct_mat)
      setProgress(1)
    })
    
    res <- run_normfinder_logic(ct_mat)
    # Abdeckung hinzufügen
    coverage <- gene_presence %>% 
      mutate(Coverage = paste0(round((n_plates / num_plates) * 100), "%")) %>%
      select(Gene, Coverage)
    
    res <- res %>% left_join(coverage, by = "Gene")
    
    return(res)
  })
  

  
  output$ref_stability_plot_multi <- renderPlotly({
    res <- validationResultsMulti()
    req(res)
    
    # Wir erstellen den Plot analog zum Single-Plate Plot, fügen aber die Coverage-Information hinzu.
    p <- ggplot(res, aes(x = reorder(Gene, Stability), y = Stability, fill = Stability,
                         text = paste0("Gene: ", Gene, 
                                       "<br>Stability: ", round(Stability, 4),
                                       "<br>Plate Coverage: ", Coverage))) +
      geom_bar(stat = "identity", alpha = 0.8) +
      # Farbskala: Grün (stabil) zu Rot (instabil)
      scale_fill_gradient(low = "#2ecc71", high = "#e74c3c") +
      theme_minimal() +
      labs(x = "Candidate Genes", 
           y = "Stability Value (lower is better)",
           title = "Multi-Plate NormFinder Ranking") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    # Tooltip auf unseren custom 'text' im aes() lenken
    ggplotly(p, tooltip = "text")
  })

  output$ref_ranking_table_multi <- renderDT({
    res <- validationResultsMulti()
    req(res)
    
    datatable(res, rownames = FALSE) %>%
      formatRound("Stability", 3) %>%
      formatStyle(
        'Coverage',
        color = styleEqual("100%", "green", "red"),
        fontWeight = styleEqual("100%", "normal", "bold")
      )
  })
  
  
  
  
  # --- NORMFINDER LOGIK ---
  
  # Hilfsfunktion zur Berechnung (basiert auf dem Andersen-Algorithmus)
  run_normfinder_logic <- function(ct_matrix, group_vector = NULL) {
    # ct_matrix: Zeilen = Proben, Spalten = Gene
    # group_vector: Optional für Inter-Group Variabilität
    
    if (is.null(group_vector)) group_vector <- rep(1, nrow(ct_matrix))
    

    genes <- colnames(ct_matrix)
    stability_scores <- sapply(genes, function(g) {
      # Berechne die Abweichung des Gens zum Mittelwert aller anderen Kandidaten
      gene_data <- ct_matrix[, g]
      other_genes <- ct_matrix[, colnames(ct_matrix) != g]
      if (is.vector(other_genes)) {
        ref_line <- other_genes
      } else {
        ref_line <- rowMeans(other_genes, na.rm = TRUE)
      }
      
      # Differenz berechnen
      diff <- gene_data - ref_line
      # Stabilität als SD der Differenzen
      sd(diff, na.rm = TRUE)
    })
    
    res <- data.frame(
      Gene = genes,
      Stability = as.numeric(stability_scores)
    ) %>% arrange(Stability) %>%
      mutate(Rank = row_number())
    
    return(res)
  }  
  

  #Reference Gene Finder Apply Button Single Plate
  observeEvent(input$apply_best_ref, {
    res <- validationResults()
    req(res)
    
    # Nimm die Top 2 Gene
    best_genes <- res$Gene[1:2]
    
    # Update das Haupt-Input für Referenzgene (aus dem Analyse-Tab)
    updateVirtualSelect("Refs", selected = best_genes)
    
    showNotification(paste("Applied", paste(best_genes, collapse = ", "), "as reference genes."), type = "message")
  }) 
  
  #Reference Gene Finder Apply Button Multi Plate
  observeEvent(input$apply_best_ref_multi, {
    res <- validationResultsMulti()
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
  
  
  
  output$algo_description <- renderUI({
    tagList(
      tags$h5("NormFinder Algorithm", style = "color: #2c3e50;"),
      tags$p("NormFinder is a model-based approach to identify the optimal normalization gene(s) among a set of candidate genes."),
      
      tags$div(
        style = "background-color: #f8f9fa; border-left: 5px solid #2ecc71; padding: 15px; margin: 10px 0;",
        tags$b("Interpretation:"),
        tags$ul(
          tags$li(tags$b("Stability Value:"), " Represents the estimated intra- and inter-group variation. "),
          tags$li(tags$b("Rule of Thumb:"), " A lower value indicates higher expression stability."),
          tags$li("Values ", tags$b("< 0.15"), " are generally considered excellent for stable reference genes.")
        )
      ),
      
      tags$h6("How it works:"),
      tags$p("Unlike pairwise methods, NormFinder evaluates each gene's stability by comparing its variation to the average of the other genes. It is particularly robust because it can account for experimental groups (e.g., Treatment vs. Control), ensuring that a gene isn't just stable in one group but across the entire experimental setup."),
      
      tags$hr(),
      
      tags$div(
        style = "font-size: 0.85rem; color: #7f8c8d;",
        tags$p(tags$i("Reference: Andersen CL, Jensen JL, Orntoft TF. Normalization of real-time quantitative reverse transcription-PCR data: a model-based variance estimation approach. Cancer Res. 2004."))
      )
    )
  })
  
  output$algo_description_multi <- renderUI({
    tagList(
      tags$h5("NormFinder Algorithm", style = "color: #2c3e50;"),
      tags$p("NormFinder is a model-based approach to identify the optimal normalization gene(s) among a set of candidate genes."),
      
      tags$div(
        style = "background-color: #f8f9fa; border-left: 5px solid #2ecc71; padding: 15px; margin: 10px 0;",
        tags$b("Interpretation:"),
        tags$ul(
          tags$li(tags$b("Stability Value:"), " Represents the estimated intra- and inter-group variation. "),
          tags$li(tags$b("Rule of Thumb:"), " A lower value indicates higher expression stability."),
          tags$li("Values ", tags$b("< 0.15"), " are generally considered excellent for stable reference genes.")
        )
      ),
      
      tags$h6("How it works:"),
      tags$p("Unlike pairwise methods, NormFinder evaluates each gene's stability by comparing its variation to the average of the other genes. It is particularly robust because it can account for experimental groups (e.g., Treatment vs. Control), ensuring that a gene isn't just stable in one group but across the entire experimental setup."),
      
      tags$hr(),
      
      tags$div(
        style = "font-size: 0.85rem; color: #7f8c8d;",
        tags$p(tags$i("Reference: Andersen CL, Jensen JL, Orntoft TF. Normalization of real-time quantitative reverse transcription-PCR data: a model-based variance estimation approach. Cancer Res. 2004."))
      )
    )
  })
  
  
  
  
  
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
