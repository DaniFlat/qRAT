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

server <- function(input, output, session) {

  # fully shutdown qRAT and R when user session ends
  session$onSessionEnded(stopApp)

  # loading screen to not have the screen flash bright white
  hostess <- Hostess$new("loader")

  for (i in 1:10) {
    Sys.sleep(runif(1) / 2)
    hostess$set(i * 10)
  }
  waiter_hide()


  # Button UpdateCheck
  observeEvent(input$updateCheck, {

    show_alert(
              title = "Check for Updates",
              text = tags$span("Visit the website to check for updates.", tags$br(), tags$br(), tags$a(href = "https://www.uibk.ac.at/de/microbiology/services/qrat/", "Go to the website")),
              html = TRUE,
              type = "info",
              showCloseButton = TRUE
            )
  })

  # Button Citation
  observeEvent(input$Citation, {
    showModal(modalDialog(
      title = "Citation",
      tags$span("If you use qRAT in published research, please cite the following paper:",
      tags$a(href = "https://doi.org/10.1186/s12859-022-04823-7", "qRAT"), tags$br(),
      "If you use ddCq values, please also cite the paper for:",
      tags$a(href = "https://doi.org/doi:10.18129/B9.bioc.ddCt", "ddCt"), tags$br(),
      "If you use the statistical analysis, please also cite the paper for:",
      tags$a(href = "https://doi.org/doi:10.18129/B9.bioc.limma", "limma"), tags$br()
      ),
      footer = tagList(
        modalButton("Close (Esc)"))
    ))
  })

  # Button housekeepingInfo
  observeEvent(input$housekeepingInfo, {
    showModal(modalDialog(
      title = "Reference Genes",
      tags$span("The stability of all appointed reference genes needs to be validated in advance. Popular algorithms to determine the most stable reference (housekeeping) genes
			from a set of candidate reference are geNorm, BestKeeper and NormFinder. On average 2-4 reference genes should ideally be used for final normalization in a given experiment.",
        tags$br(),
        tags$a(href = "https://doi.org/10.1186%2Fgb-2002-3-7-research0034", "More Information")
      ),
      footer = tagList(
        modalButton("Close (Esc)"))
    ))
  })

  # Button housekeepingInfo Multiple Plates
  observeEvent(input$housekeepingInfoMP, {
    showModal(modalDialog(
      title = "Reference Genes",
      tags$span("The stability of all appointed reference genes needs to be validated in advance. Popular algorithms to determine the most stable reference (housekeeping) genes
			from a set of candidate reference are geNorm, BestKeeper and NormFinder. On average 2-4 reference genes should ideally be used for final normalization in a given experiment.",
                tags$br(),
                tags$a(href = "https://doi.org/10.1186%2Fgb-2002-3-7-research0034", "More Information")
      ),
      footer = tagList(
        modalButton("Close (Esc)"))
    ))
  })

  # Button outlierSP
  observeEvent(input$outlierSP, {
    showModal(modalDialog(
      title = "Define outlier",
      tags$span("WIP define outlier; Wells marked as outlier will be neglected in subsequent filtering and analysis",
                tags$br()
      ),
      DTOutput("defineOutlierSP"),
      footer = tagList(
        modalButton("Close (Esc)"))
    ))
  })


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

  SP_validateInputFile <- function(Sample_Column, Well_Column, Gene_Column, Cq_Column) {
    modalDialog(
      title = "Validate Single Plate File",

        div(icon("circle-check", style = "color: green;"), tags$b("File type is correct.", style = "color: green;")),

      if (Sample_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Sample Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Sample Column is missing...", style = "color: red;"))
      },
      if (Well_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Well Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Well Column is missing...", style = "color: red;"))
      },
      if (Gene_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Gene Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Gene Column is missing...", style = "color: red;"))
      },
      if (Cq_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Cq Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Cq Column is missing...", style = "color: red;"))
      },

      footer = tagList(
        modalButton("Close (Esc)")
      )
    )
  }

  MP_validateInputFile <- function(Sample_Column, Well_Column, Gene_Column, Cq_Column, plateName) {
    modalDialog(
      title = "Validate Multiple Plate Files",

    #  textOutput("plateName"), #this does not work yet; should show the filename of the respective plate

      div(icon("circle-exclamation"), tags$b("Files will be checked one by one.")),

      div(icon("circle-check", style = "color: green;"), tags$b("File type is correct.", style = "color: green;")),

      if (Sample_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Sample Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Sample Column is missing...", style = "color: red;"))
      },
      if (Well_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Well Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Well Column is missing...", style = "color: red;"))
      },
      if (Gene_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Gene Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Gene Column is missing...", style = "color: red;"))
      },
      if (Cq_Column){
        div(icon("circle-check", style = "color: green;"), tags$b("Found Cq Column", style = "color: green;"))}
      else {
        div(icon("circle-exclamation", style = "color: red;"), tags$b("Cq Column is missing...", style = "color: red;"))
      },

      footer = tagList(
        modalButton("Close (Esc)")
      )
    )
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



  # Single Plate Data Manipulations (filtering, quality, preparations for subsequent calculations)

  setData <- reactive({
    req(input$dtfile)
    info <- readSinglePlateData()
    xx <- info$data
    if (is.null(xx)) {
      return(NULL)
    }

    if (input$Sep == "tab") sep <- "\t"
    if (input$Sep == "comma") sep <- ","
    if (input$Sep == "semicolon") sep <- ";"
    if (input$Sep == "auto") sep <- "\t"

    experimental_controls <- input$NTC_Input




    #####manual definition of outliers

    # create a character vector of shiny inputs
    shinyInput = function(FUN, len, id, value, ...) {
      if (length(value) == 1) value <- rep(value, len)
      inputs = character(len)
      for (i in seq_len(len)) {
        inputs[i] = as.character(FUN(paste0(id, i), label = NULL, value = value[i]))
      }
      inputs
    }

    # obtain the values of inputs
    shinyValue = function(id, len) {
      unlist(lapply(seq_len(len), function(i) {
        value = input[[paste0(id, i)]]
        if (is.null(value)) TRUE else value
      }))
    }

    #prepare data for outlier definition
    markOutlierTable <- xx
    n <- nrow(markOutlierTable)
    markOutlierTable %>%
      mutate(outlier = shinyInput(checkboxInput, n, 'cb_', value = TRUE, width='1px')) %>%
      mutate(YN = rep(TRUE, n))

    ## Outlier Definition SP
    loopData = reactive({
      markOutlierTable$cb <<- shinyInput(checkboxInput, n, 'cb_', value = shinyValue('cb_', n), width='1px')
      markOutlierTable$YN <<- shinyValue('cb_', n)
      markOutlierTable
    })

    output$defineOutlierSP = DT::renderDataTable(
      isolate(loopData()),
      escape = FALSE, selection = 'none',
      options = list(
        dom = 't', paging = FALSE, ordering = FALSE,
        preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); } ')
      ))


    # remove (set to NA) replicates based on threshold (bad reps and max ct) and empty rows (no name, no gene)
    correctRep <- xx %>%
      group_by(Sample, Gene) %>%
      mutate(Cq_backup = Ct) %>%
      mutate(Ct = remove_replicates(Ct)) %>%
      mutate(Ct = remove_maxCt(Ct)) %>%
      filter(Sample != "") %>%
      filter(!Sample %in% experimental_controls) %>%
      filter(Gene != "") %>%
      ungroup()



    badRep <- correctRep %>%
      filter(is.na(Ct)) %>%
      rename(Cq = Cq_backup)

    xx <- correctRep


    fwrite(xx, file = "fileSingle.csv", sep = sep)

    Genes <- unique(xx$Gene)
    refs <- grep("^(ACT|UBQ|REF|RF|BACTIN|ACTIN|GAPDH|TUB|TUBULIN).*", Genes, ignore.case = TRUE)
    if (length(refs) == 0) {
      refs <- Genes[length(Genes)]
    } else {
      refs <- Genes[refs]
    }

    # Update Inputfield for Gene Input
    updateVirtualSelect("Refs", choices = Genes, selected = refs)

    ## read to qPCR format single plate
    htset <- NULL
    ddset <- NULL
    xdata <- read.qPCRtable("fileSingle.csv", sep = sep)
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))

      # Update Inputfields (Plot Settings and Input Settings)

      updateVirtualSelect("Mock", choices = Samples, selected = Samples[1])
      updateVirtualSelect("Comps", choices = Samples, selected = Samples) # update input statistics

      if (length(Samples) > 4 | max(nchar(Samples)) > 6) {
        updateSliderInput(session, "side1", value = 8)
        updateSliderInput(session, "xsrt", value = 30)
      }
    }

    return(list(data = xx, data.ht = htset, data.ddct = ddset, badReplicates = badRep, CompSamples = Samples))
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



  # read Multiple Plate Data

  readMPData <- reactive({
    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"
    if (input$SepM == "auto") sep <- "auto"

    req(input$plates)

    f <- input$plates$datapath
    ff <- input$plates
    if (is.null(f)) {
      return(NULL)
    } else {

      # access single files (get number of files before) and check header
      numfiles <- nrow(ff)
    }
    filelist <- list()

    #Check file type
    req(input$plates)
    for (i in 1:numfiles) {
      ext <- file_ext(ff[[i, "datapath"]])
    if (!ext %in% c(
      'text/csv',
      'text/comma-separated-values',
      'text/tab-separated-values',
      'text/plain',
      'csv',
      'txt',
      'tsv')) {
      showModal(check_filetype(failed = TRUE))
    }}

    for (i in 1:numfiles) {
      file <- ff[[i, "datapath"]]
      file.header <- readLines(con = file, n = 100)
      n.header <- grep("^Well", file.header) - 1
      if (length(n.header) == 0) {
        n.header <- 0
      }
      filelist[[i]] <- fread(file, header = TRUE, sep = sep, skip = n.header)


      ##Validate Input Data File
      #change all column names to lowercase for easy comparison with prepared column-name-lists
      column_namesTest <- tolower(names(filelist[[i]]))

      validateInputFile <- function(Sample_Column=TRUE, Well_Column=TRUE, Gene_Column=TRUE, Cq_Column=TRUE, plateName) {
        if (length(intersect(sampleColumnList,column_namesTest)) < 1) Sample_Column = FALSE
        if (length(intersect(wellColumnList,column_namesTest)) < 1) Well_Column = FALSE
        if (length(intersect(geneColumnList,column_namesTest)) < 1) Gene_Column = FALSE
        if (length(intersect(cqColumnList,column_namesTest)) < 1) Cq_Column = FALSE

        plateName <- ff$name

        showModal(MP_validateInputFile(Sample_Column, Well_Column, Gene_Column, Cq_Column, plateName))
      }

      validateInputFile()
    }

    multiplePlates <- rbindlist(filelist, use.names = TRUE, fill = TRUE, idcol = "PlateNumber")

    #change all column names to lowercase for easy comparison with prepared column-name-lists
    multiplePlates <- multiplePlates %>%
      rename_with(tolower)

    ## change all column names to lowercase and rename columns to Well, Sample, Gene, Ct
    multiplePlates <- multiplePlates %>%
      rename(
        Sample = any_of(sampleColumnList),
        Well = any_of(wellColumnList),
        Gene = any_of(geneColumnList),
        Ct = any_of(cqColumnList),
        PlateNumber = platenumber
      )


    ## generate replicate number
    multiplePlates$TempRepNum <- paste(multiplePlates$Sample, multiplePlates$Gene)
    multiplePlates$rp.num <- ave(multiplePlates$Sample, multiplePlates$TempRepNum, FUN = seq_along)
    multiplePlates <- subset(multiplePlates, select = c(PlateNumber, Well, Sample, Gene, Ct, rp.num)) %>%
      filter(Well != "")


    ## replace "," with "."
    CorrectCt <- multiplePlates$Ct
    CorrectCt <- gsub(",", ".", CorrectCt, ignore.case = TRUE)
    multiplePlates$Ct <- as.numeric(CorrectCt)

    # Get Plate Number to Input Selection for Spatial Plate View
    ChosenPlate <- unique(as.character(multiplePlates$PlateNumber))
    updateVirtualSelect("PlateSelect", choices = ChosenPlate, selected = ChosenPlate[1])

    # raw data for plots (MP spatial plate view; MP cq distribution)
    MultiplePlatesView <- multiplePlates %>%
      separate(Well,
               into = c("text", "num"),
               sep = "(?<=[A-Za-z])(?=[0-9])"
      )

    # Update Inputfield for MP_NTC Selection with all sample names
    MP_NTCs <- unique(as.character(multiplePlates$Sample))
    updateVirtualSelect("NTC_Input_MP", choices = MP_NTCs)

    multiplePlates

    return(list(data = multiplePlates, MPV = MultiplePlatesView))
  })

  # Multiple Plate read files

  multiData <- reactive({

    req(input$plates)
    info <- readMPData()
    multiplePlates <- info$data
    if (is.null(multiplePlates)) {
      return(NULL)
    }

    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"
    if (input$SepM == "auto") sep <- "\t"

    experimental_controls <- input$NTC_Input_MP

    ## remove (set to NA) replicates based on threshold (bad reps and max ct) and empty rows (no name, no gene)
    correctRep <- multiplePlates %>%
      group_by(PlateNumber, Sample, Gene) %>%
      mutate(Cq_backup = Ct) %>%
      mutate(Ct = remove_replicatesMP(Ct)) %>%
      mutate(Ct = remove_maxCtMulti(Ct)) %>%
      filter(Sample != "") %>%
      filter(!Sample %in% experimental_controls) %>%
      filter(Gene != "") %>%
      ungroup()

    badRep <- correctRep %>%
      filter(is.na(Ct)) %>%
      rename(Cq = Cq_backup)

    multiplePlates <- correctRep

    ## read Samples for IPC extraction
    multiplePlatesCOPY <- multiplePlates
    SamplesWithIPCs <- unique(as.character(multiplePlatesCOPY$Sample))
    updateVirtualSelect("IPC", choices = SamplesWithIPCs, selected = SamplesWithIPCs[0])


    fwrite(multiplePlates, file = "fileMulti.csv", sep = sep)


    return(list(data = multiplePlates, badReplicates = badRep))
  })


  ## multiData 2, for separating IPCs from main Samples

  multiData2 <- reactive({
    req(input$plates)
    # req(input$IPC)

    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"
    if (input$SepM == "auto") sep <- "\t"

    info <- multiData()
    if (is.null(info)) {
      return(NULL)
    }
    info2 <- info$data
    if (is.null(info2)) {
      return(NULL)
    }
    if (input$Id027 == TRUE) {
      req(input$IPC)
      ipcs <- input$IPC

      # remove IPCs from Samples, as they interfere with expression calculations
      info2 <- info2[!(info2$Sample == ipcs), ]

      fwrite(info2, file = "fileMultiNoIPCs.csv", sep = sep)
    } else {

      # don't try to remove IPCs from samples
      fwrite(info2, file = "fileMultiNoIPCs.csv", sep = sep)
    }


    ## read genes
    Genes <- unique(info2$Gene)
    refs <- grep("^(ACT|UBQ|REF|RF|BACTIN|ACTIN|GAPDH|TUB|TUBULIN).*", Genes, ignore.case = TRUE)
    if (length(refs) == 0) {
      refs <- Genes[length(Genes)]
    } else {
      refs <- Genes[refs]
    }
    updateVirtualSelect("RefsM", choices = Genes, selected = refs)



    ## read to qPCR format multiple plates
    htset <- NULL
    ddset <- NULL
    xdata <- read.qPCRtableMulti("fileMultiNoIPCs.csv", sep = sep)
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))

      # Update Inputfields (Sample Selection)

      updateVirtualSelect("MockM", choices = Samples, selected = Samples[1])
      updateVirtualSelect("CompsM", choices = Samples, selected = Samples) # update input statistics

      if (length(Samples) > 4 | max(nchar(Samples)) > 6) {
        updateSliderInput(session, "side1", value = 8)
        updateSliderInput(session, "xsrt", value = 30)
      }
    }

    return(list(SamplesNoIPCs = info2, xdata, data.ht = htset, data.ddct = ddset, CompSamples = Samples))
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
    updateOrderInput(session, "x_order_MultidCq", items = SamplesdCtMulti)
    updateVirtualSelect("SamplePickerDDCtMulti", choices = SamplesddCtMulti, selected = SamplesddCtMulti)
    updateOrderInput(session, "x_order_MultiddCq", items = SamplesddCtMulti)
    updateVirtualSelect("GenePickerMulti", choices = GenesdCtMulti, selected = GenesdCtMulti)
    updateVirtualSelect("GenePickerDDCtMulti", choices = GenesddCtMulti, selected = GenesddCtMulti)

    ## updateRadioButtons(session, 'geneNameRel', choices=targets, selected=targets[1])
    ## updateCheckboxGroupInput(session, 'geneNameAbs', choices=targets, selected=targets)
    list(multiRelative = expr.rel, multiAbsolute = expr.abs, elist = resultElist)
  })




  ## extract IPCs; IPCs must have the same Sample name across plates

  extractIPCs <- reactive({
    if (input$Id027 == TRUE) {
      req(input$IPC)

      info <- multiData()
      if (is.null(info)) {
        return(NULL)
      }
      info2 <- info$data
      if (is.null(info2)) {
        return(NULL)
      }
      ipcs <- input$IPC
      if (is.null(ipcs)) {
        return(NULL)
      }

      if (input$IPC == "") {
        return(NULL)
      }

      extractedIPCs <- info2[info2$Sample %in% ipcs, ]

      extractedIPCs$Sample <- paste(extractedIPCs$Sample, extractedIPCs$PlateNumber)

      # generate new replicate number
      extractedIPCs$TempRepNum <- paste(extractedIPCs$Sample, extractedIPCs$Gene)
      extractedIPCs$rp.num <- ave(extractedIPCs$Sample, extractedIPCs$TempRepNum, FUN = seq_along)

      # normalize IPCs

      if (input$SepM == "tab") sep <- "\t"
      if (input$SepM == "comma") sep <- ","
      if (input$SepM == "semicolon") sep <- ";"
      if (input$SepM == "auto") sep <- "\t"

      extractedIPCs <- na.omit(extractedIPCs)
      fwrite(extractedIPCs, file = "IPCs.csv", sep = sep)

      gm_mean <- function(x, na.rm = TRUE) {
        exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
      }

      geometricMean <- gm_mean(extractedIPCs$Ct)

      extractedIPCs <- extractedIPCs %>%
        mutate(CtDivGmean = Ct / geometricMean)

      CalibrationFactors <- extractedIPCs %>%
        group_by(PlateNumber) %>%
        summarize(Mean = mean(CtDivGmean, na.rm = TRUE))


      return(list(extrIPCs = extractedIPCs, CalibrationFactors = CalibrationFactors))
    } else {

    }
  })



  # calibrate Ct values with calibration factor, inter-plate calibration
  interPlateCalibration <- reactive({
    if (input$Id027 == TRUE) {
      calF <- extractIPCs()
      if (is.null(calF)) {
        return(NULL)
      }
      data <- multiData2()
      if (is.null(data)) {
        return(NULL)
      }

      x <- data$SamplesNoIPCs

      x %>%
        left_join(calF$CalibrationFactors, by = "PlateNumber") %>%
        mutate(Ct = Ct / Mean) %>%
        select(-Mean)
    } else {
      data <- multiData2()
      if (is.null(data)) {
        return(NULL)
      }
      x <- data$SamplesNoIPCs
    }
  })


  interPlateCalibration_Plot <- reactive({

    info <- interPlateCalibration()
    if (is.null(info)) NULL

    info2 <- multiData2()
    if (is.null(info2)) NULL

      calibrated <- info %>% unite("Sample_Gene", Sample:Gene, remove = FALSE)
      noncalibrated <- info2$SamplesNoIPCs %>% unite("Sample_Gene", Sample:Gene, remove = FALSE)

      SamplesIPCcomparison <- unique(calibrated$Sample_Gene)
      updateVirtualSelect("SamplePickerIPCcomparison", choices = SamplesIPCcomparison, selected = "")

      return(list(calibrated = calibrated, noncalibrated = noncalibrated))


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



  ## dCt Expression Graph Single Plate

  output$ddctAbsGraph <- renderPlotly({
    colorPick <- input$colorpicker
    scalePick <- input$scale
    formatChoice <- input$exportFormat
    PlotDataPick <- input$PlotData
    PlotType <- input$PlotType
    PlotWidth <- input$width_dCq
    PlotHeight <- input$height_dCq
    PlotScale <- input$scale_dCq
    plotTitle <- input$plotTitle_dCq
    legendTitle <- input$legendTitle_dCq
    legendPosition <- input$legendPosition_dCq
    xTextAngle <- input$xTextAngle_dCq

    info <- resultDdct()
    if (is.null(info)) {
      return(NULL)
    }

    # get the correct data for plot depending on user input (RQ or dCq values)
    if (PlotDataPick == "RQ") {
      PlotDataError <- "expr.sd"
      PlotDataPick <- "expr"
      yTitle <-"relative quantity"
    } else {
      if (PlotDataPick == "-dCq") {
        PlotDataError <- "dCt.sd"
        PlotDataPick <- "negdCq"
        yTitle <- "-\u0394Cq"
      } else {
        PlotDataError <- "dCt.sd"
        PlotDataPick <- "dCt"
        yTitle <- "\u0394Cq"
      }
    }

    df <- info$absolute %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% input$SamplePicker) %>%
      filter(Gene %in% input$GenePicker)

    df <- rename(df, "negdCq" = "-dCq")


    # Convert x-axis variable to factor with desired levels and order
    df$Sample <- factor(df$Sample, levels = input$x_order)

    # Define color palette based on user input and unique Gene levels
    numberOfGenes <- n_distinct(df$Gene)
    GeneNames <- pull(distinct(df, Gene))
    GeneNames <- droplevels(GeneNames)


    if (colorPick == "viridis") {
      colourpalette <- viridis(numberOfGenes)
      names(colourpalette) <- GeneNames
    } else {
      colourpalette <- colorRampPalette(brewer.pal(numberOfGenes, colorPick))(numberOfGenes)
      names(colourpalette) <- GeneNames

    }


  if (PlotType == "Bar Chart") {
   figNormal <- ggbarplot(df, x="Sample",y=PlotDataPick,
              fill = "Gene", color = "Gene", palette = colourpalette,
              position = position_dodge2(0.9, preserve = "single")) + aes(text = paste("Sample:", Sample,
                                                                 "<br>Gene:", Gene,
                                                                 "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                 "<br>Error:", round(.data[[PlotDataError]], 3)))


    } else {

      figNormal <- ggdotchart(df, x="Sample",y=PlotDataPick,
                             fill = "Gene", color = "Gene", palette = colourpalette, sorting = "none",
                             position = position_dodge(0.9), size=3) + aes(text = paste("Sample:", Sample,
                                                                                        "<br>Gene:", Gene,
                                                                                        "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                        "<br>Error:", round(.data[[PlotDataError]], 3))) + guides(fill = guide_legend(override.aes = list(color = NULL)))

    }


    if(input$ShowErrorBar==TRUE){
      if(PlotType == "Bar Chart") {

        figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                   ymax = get(PlotDataPick) + get(PlotDataError)),
                                               width=.7, color = "black", position = position_dodge2(0.9),
                                               show.legend = FALSE)
      } else {
        figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                   ymax = get(PlotDataPick) + get(PlotDataError)),
                                               width=.7, color = "black", position = position_dodge(0.9),
                                               show.legend = FALSE)
          }
      }

    #Additional Plot Settings ggpubr
    figNormal <- ggpar(figNormal, x.text.angle = xTextAngle, legend = legendPosition, title = plotTitle, legend.title = legendTitle)

    if (scalePick == "normal") {
        fig <- ggplotly(figNormal, tooltip = c("text")) %>%
          style(hoverlabel = list(bgcolor = "white"))
      } else {
        fig <- figNormal + yscale("log2", .format = TRUE)
        fig <- ggplotly(fig, tooltip = c("text")) %>%
          style(hoverlabel = list(bgcolor = "white"))
      }

    fig <- fig %>% layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t)

    #Plotly Settings and Plot Export
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = FALSE, displayModeBar = TRUE
      )

    fig
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

  ## ddCt Expression Graph Single Plate

  output$ddctRelGraph <- renderPlotly({
    info <- resultDdct()
    if (is.null(info)) {
      return(NULL)
    }
    colorPick <- input$colorpickerDDCt
    scalePick <- input$scaleDDCt
    formatChoice <- input$exportFormatDDCt
    PlotDataPick <- input$PlotDataDDCt
    PlotType <- input$PlotTypeDDCt
    PlotWidth <- input$width_ddCq
    PlotHeight <- input$height_ddCq
    PlotScale <- input$scale_ddCq
    plotTitle <- input$plotTitle_ddCq
    legendTitle <- input$legendTitle_ddCq
    legendPosition <- input$legendPosition_ddCq
    xTextAngle <- input$xTextAngle_ddCq

    SamplePicker <- input$SamplePickerDDCt
    GenePicker <- input$GenePickerDDCt

    df <- info$relative %>%
        mutate(Sample = as.character(Sample)) %>%
        filter(Sample %in% c(SamplePicker)) %>%
        filter(Gene %in% c(GenePicker))

    # get the correct data for plot depending on user input (FC or ddCT values)
    if (PlotDataPick == "Fold Change") {
      PlotDataError <- "FC.sd"
      PlotDataPick <- "FC"
      yTitle <- "fold change"
    } else {
      PlotDataError <- "ddCt.sd"
      PlotDataPick <- "ddCt"
      yTitle <- "\u0394\u0394Cq"
    }

    # Convert x-axis variable to factor with desired levels and order
    df$Sample <- factor(df$Sample, levels = input$x_orderDDCt)

    # Define color palette based on user input and unique Gene levels
    numberOfGenes <- n_distinct(df$Gene)
    GeneNames <- pull(distinct(df, Gene))
    GeneNames <- droplevels(GeneNames)

    if (colorPick == "viridis") {
      colourpalette <- viridis(numberOfGenes)
      names(colourpalette) <- GeneNames
    } else {
      colourpalette <- colorRampPalette(brewer.pal(numberOfGenes, colorPick))(numberOfGenes)
      names(colourpalette) <- GeneNames

    }


    if (PlotType == "Bar Chart") {
      figNormal <- ggbarplot(df, x="Sample",y=PlotDataPick,
                             fill = "Gene", color = "Gene", palette = colourpalette,
                             position = position_dodge2(0.9, preserve = "single")) + aes(text = paste("Sample:", Sample,
                                                                                "<br>Gene:", Gene,
                                                                                "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                "<br>Error:", round(.data[[PlotDataError]], 3)))


    } else {

      figNormal <- ggdotchart(df, x="Sample",y=PlotDataPick,
                              fill = "Gene", color = "Gene", palette = colourpalette, sorting = "none",
                              position = position_dodge(0.9), size=3) + aes(text = paste("Sample:", Sample,
                                                                                         "<br>Gene:", Gene,
                                                                                         "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                         "<br>Error:", round(.data[[PlotDataError]], 3))) + guides(fill = guide_legend(override.aes = list(color = NULL)))

    }


    if(input$ShowErrorBarDDCt==TRUE){
      if (PlotType == "Bar Chart") {
          figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                 ymax = get(PlotDataPick) + get(PlotDataError)),
                                             width=.7, color = "black", position = position_dodge2(0.9),
                                             show.legend = FALSE)
      } else {
        figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                   ymax = get(PlotDataPick) + get(PlotDataError)),
                                               width=.7, color = "black", position = position_dodge(0.9),
                                               show.legend = FALSE)
      }
    }

    #Additional Plot Settings ggpubr
    figNormal <- ggpar(figNormal, x.text.angle = xTextAngle, legend = legendPosition, title = plotTitle, legend.title = legendTitle)

    if (scalePick == "normal") {
      fig <- ggplotly(figNormal, tooltip = c("text")) %>%
        style(hoverlabel = list(bgcolor = "white"))
    } else {
      fig <- figNormal + yscale("log2", .format = TRUE)
      fig <- ggplotly(fig, tooltip = c("text")) %>%
        style(hoverlabel = list(bgcolor = "white"))
    }

    fig <- fig %>% layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t)

    #Plotly Settings and Plot Export
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = FALSE, displayModeBar = TRUE
      )

    fig
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


  ## Comparison of calibrated data with non-calibrated data

  output$comparisonCalibrated <- renderPlotly({
    if (input$Id027 == TRUE) {
      info <- interPlateCalibration_Plot()
      if (is.null(info)) NULL

      calibrated <- info$calibrated  %>%
        group_by(Sample_Gene, PlateNumber) %>%
        summarize(MeanCq = mean(Ct, na.rm=TRUE))
      calibrated_sd <- info$calibrated  %>%
        group_by(Sample_Gene, PlateNumber) %>%
        summarize(sdCq = sd(Ct, na.rm = TRUE))
      calibrated <- data.frame(calibrated, calibrated_sd$sdCq)
      calibrated <- rename(calibrated, c("sdCq" = "calibrated_sd.sdCq"))


      noncalibrated <- info$noncalibrated  %>%
        group_by(Sample_Gene, PlateNumber) %>%
        summarize(MeanCq = mean(Ct, na.rm=TRUE))
      noncalibrated_sd <- info$noncalibrated  %>%
        group_by(Sample_Gene, PlateNumber) %>%
        summarize(sdCq = sd(Ct, na.rm = TRUE))
      noncalibrated <- data.frame(noncalibrated, noncalibrated_sd$sdCq)
      noncalibrated <- rename(noncalibrated, c("sdCq" = "noncalibrated_sd.sdCq"))




      #SamplePicker for plotting the comparison
      SamplePicker <- input$SamplePickerIPCcomparison

      #create additional column "label" for labeling/coloring the bars in subsequent plot correctly
      calibrated$label <- NA
      calibrated$label <- paste("Plate", calibrated$PlateNumber, "calibrated")
      noncalibrated$label <- NA
      noncalibrated$label <- paste("Plate", noncalibrated$PlateNumber, "non-calibrated")

      calibrated <- calibrated %>% filter(Sample_Gene %in% c(SamplePicker))
      noncalibrated <- noncalibrated %>% filter(Sample_Gene %in% c(SamplePicker))

      fig1 <- plot_ly(colors = "Spectral", marker = list(size = 10, line = list(color = "rgba(0, 0, 0, .8)", width = 2))) %>%
        add_trace(data = as.data.frame(calibrated),
        x = ~Sample_Gene, y = ~MeanCq, color = ~as.character(label), type = "bar") %>%
        add_trace(data = as.data.frame(noncalibrated), x = ~Sample_Gene, y = ~MeanCq, color = ~as.character(label), type ="bar") %>%
        layout(
          title = "Comparison calibrated/non-calibrated",
          font=t,
          xaxis = xaxis,
          yaxis = yaxis,
          barmode = 'group',
          bargroupgap = 0.1,
          margin = list(t = 50),
          legend=list(title=list(text='Plates'))
        ) %>%
        config(
          displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
          toImageButtonOptions = list(format = "png", scale = 3, width=768, height=768), scrollZoom = TRUE, displayModeBar = TRUE
        )  %>%
        layout(yaxis = list(title = "Mean Cq"), xaxis = list(title = "Sample"), font=t)

    } else {}
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

  ## dCt Gene Expression Graph Multiple Plates

  output$ddctMultiAbsGraph <- renderPlotly({
    info <- resultDdctMulti()
    if (is.null(info)) {
      return(NULL)
    }

    colorPick <- input$colorpickerMulti
    scalePick <- input$scaleMulti
    formatChoice <- input$exportFormatMulti
    SamplePicker <- input$SamplePickerMulti
    GenePicker <- input$GenePickerMulti
    PlotDataPick <- input$PlotDataMulti
    PlotType <- input$PlotTypeMulti
    PlotWidth <- input$width_dCqMulti
    PlotHeight <- input$height_dCqMulti
    PlotScale <- input$scale_dCqMulti
    plotTitle <- input$plotTitle_MultidCq
    legendTitle <- input$legendTitle_MultidCq
    legendPosition <- input$legendPosition_MultidCq
    xTextAngle <- input$xTextAngle_MultidCq


    # get the correct data for plot depending on user input (RQ or dCT values)
    if (PlotDataPick == "RQ") {
      PlotDataError <- "expr.sd"
      PlotDataPick <- "expr"
      yTitle <- "relative quantity"

    } else {
      if (PlotDataPick == "-dCq") {
        PlotDataError <- "dCt.sd"
        PlotDataPick <- "-dCq"
        yTitle <- "-\u0394Cq"
      } else {
        PlotDataError <- "dCt.sd"
        PlotDataPick <- "dCt"
        yTitle <- "\u0394Cq"
      }
    }
    #


    df <- info$multiAbsolute %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% c(SamplePicker)) %>%
      filter(Gene %in% c(GenePicker))

    df <- rename(df, "negdCq" = "-dCq")

    # Convert x-axis variable to factor with desired levels and order
    df$Sample <- factor(df$Sample, levels = input$x_order_MultidCq)

    # Define color palette based on user input and unique Gene levels
    numberOfGenes <- n_distinct(df$Gene)
    GeneNames <- pull(distinct(df, Gene))
    GeneNames <- droplevels(GeneNames)


    if (colorPick == "viridis") {
      colourpalette <- viridis(numberOfGenes)
      names(colourpalette) <- GeneNames
    } else {
      colourpalette <- colorRampPalette(brewer.pal(numberOfGenes, colorPick))(numberOfGenes)
      names(colourpalette) <- GeneNames

    }

    if (PlotType == "Bar Chart") {
      figNormal <- ggbarplot(df, x="Sample",y=PlotDataPick,
                             fill = "Gene", color = "Gene", palette = colourpalette,
                             position = position_dodge2(0.9, preserve = "single")) + aes(text = paste("Sample:", Sample,
                                                                                "<br>Gene:", Gene,
                                                                                "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                "<br>Error:", round(.data[[PlotDataError]], 3)))


    } else {

      figNormal <- ggdotchart(df, x="Sample",y=PlotDataPick,
                              fill = "Gene", color = "Gene", palette = colourpalette, sorting = "none",
                              position = position_dodge(0.9), size=3) + aes(text = paste("Sample:", Sample,
                                                                                         "<br>Gene:", Gene,
                                                                                         "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                         "<br>Error:", round(.data[[PlotDataError]], 3))) + guides(fill = guide_legend(override.aes = list(color = NULL)))

    }


    if(input$ShowErrorBarMulti==TRUE){
      if (PlotType == "Bar Chart") {
          figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                 ymax = get(PlotDataPick) + get(PlotDataError)),
                                             width=.7, color = "black", position = position_dodge2(0.9),
                                             show.legend = FALSE)
          } else {
            figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                       ymax = get(PlotDataPick) + get(PlotDataError)),
                                                   width=.7, color = "black", position = position_dodge(0.9),
                                                   show.legend = FALSE)
          }
    }

    #Additional Plot Settings ggpubr
    figNormal <- ggpar(figNormal, x.text.angle = xTextAngle, legend = legendPosition, title = plotTitle, legend.title = legendTitle)

    if (scalePick == "normal") {
      fig <- ggplotly(figNormal, tooltip = c("text")) %>%
        style(hoverlabel = list(bgcolor = "white"))
    } else {
      fig <- figNormal + yscale("log2", .format = TRUE)
      fig <- ggplotly(fig, tooltip = c("text")) %>%
        style(hoverlabel = list(bgcolor = "white"))
    }

    fig <- fig %>% layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t)

    #Plotly Settings and Plot Export
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = FALSE, displayModeBar = TRUE
      )

    fig
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

  # ddCq Gene Expression Graph Multiple Plates

  output$ddctRelGraphMulti <- renderPlotly({
    info <- resultDdctMulti()
    if (is.null(info)) {
      return(NULL)
    }
    colorPick <- input$colorpickerDDCtMulti
    scalePick <- input$scaleDDCtMulti
    formatChoice <- input$exportFormatDDCtMulti
    SamplePicker <- input$SamplePickerDDCtMulti
    GenePicker <- input$GenePickerDDCtMulti
    PlotDataPick <- input$PlotDataDDCtMulti
    PlotType <- input$PlotTypeDDCtMulti
    PlotWidth <- input$width_ddCqMulti
    PlotHeight <- input$height_ddCqMulti
    PlotScale <- input$scale_ddCqMulti
    plotTitle <- input$plotTitle_MultiddCq
    legendTitle <- input$legendTitle_MultiddCq
    legendPosition <- input$legendPosition_MultiddCq
    xTextAngle <- input$xTextAngle_MultiddCq

    # get the correct data for plot depending on user input (FC or ddCT values)
    if (PlotDataPick == "Fold Change") {
      PlotDataError <- "FC.sd"
      PlotDataPick <- "FC"
      yTitle <- "fold change"
    } else {
      PlotDataError <- "ddCt.sd"
      PlotDataPick <- "ddCt"
      yTitle <- "\u0394\u0394Cq"
    }
    #

    df <- info$multiRelative %>%
      mutate(Sample = as.character(Sample)) %>%
      filter(Sample %in% c(SamplePicker)) %>%
      filter(Gene %in% c(GenePicker))

    # Convert x-axis variable to factor with desired levels and order
    df$Sample <- factor(df$Sample, levels = input$x_order_MultiddCq)

    # Define color palette based on user input and unique Gene levels
    numberOfGenes <- n_distinct(df$Gene)
    GeneNames <- pull(distinct(df, Gene))
    GeneNames <- droplevels(GeneNames)


    if (colorPick == "viridis") {
      colourpalette <- viridis(numberOfGenes)
      names(colourpalette) <- GeneNames
    } else {
      colourpalette <- colorRampPalette(brewer.pal(numberOfGenes, colorPick))(numberOfGenes)
      names(colourpalette) <- GeneNames

    }

    if (PlotType == "Bar Chart") {
      figNormal <- ggbarplot(df, x="Sample",y=PlotDataPick,
                             fill = "Gene", color = "Gene", palette = colourpalette,
                             position = position_dodge2(0.9, preserve = "single")) + aes(text = paste("Sample:", Sample,
                                                                                "<br>Gene:", Gene,
                                                                                "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                "<br>Error:", round(.data[[PlotDataError]], 3)))


    } else {

      figNormal <- ggdotchart(df, x="Sample",y=PlotDataPick,
                              fill = "Gene", color = "Gene", palette = colourpalette, sorting = "none",
                              position = position_dodge(0.9), size=3) + aes(text = paste("Sample:", Sample,
                                                                                         "<br>Gene:", Gene,
                                                                                         "<br>Value:", round(.data[[PlotDataPick]], 3),
                                                                                         "<br>Error:", round(.data[[PlotDataError]], 3))) + guides(fill = guide_legend(override.aes = list(color = NULL)))

    }


    if(input$ShowErrorBarDDCtMulti==TRUE){
      if (PlotType == "Bar Chart") {
          figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                 ymax = get(PlotDataPick) + get(PlotDataError)),
                                             width=.7, color = "black", position = position_dodge2(0.9),
                                             show.legend = FALSE)
      } else {
        figNormal <- figNormal + geom_errorbar(aes(ymin = get(PlotDataPick) - get(PlotDataError),
                                                   ymax = get(PlotDataPick) + get(PlotDataError)),
                                               width=.7, color = "black", position = position_dodge(0.9),
                                               show.legend = FALSE)
      }
    }

    #Additional Plot Settings ggpubr
    figNormal <- ggpar(figNormal, x.text.angle = xTextAngle, legend = legendPosition, title = plotTitle, legend.title = legendTitle)

    if (scalePick == "normal") {
      fig <- ggplotly(figNormal, tooltip = c("text")) %>%
        style(hoverlabel = list(bgcolor = "white"))
    } else {
      fig <- figNormal + yscale("log2", .format = TRUE)
      fig <- ggplotly(fig, tooltip = c("text")) %>%
        style(hoverlabel = list(bgcolor = "white"))
    }

    fig <- fig %>% layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t)

    #Plotly Settings and Plot Export
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = FALSE, displayModeBar = TRUE
      )

    fig
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


}
