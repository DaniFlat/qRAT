# qPCR - Relative Expression Analysis Tool
# version: 0.1.3
#
# PLEASE CITE
# Please cite the published manuscript in all studies using qRAT
# For authors and journal information, please refer to the qRAT website https://uibk.ac.at/microbiology/services/qrat
#
#
# MIT License
#
# Copyright (c) 2022 Daniel Flatschacher
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
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
library("reshape2")
library("scales")
library("xtable")
library("data.table")
library("DT")
library("dplyr")
library("waiter")
library("thematic")
library("tidyr")
library("stringr")
library("magrittr")
library("shinycssloaders")
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
  
  thematic_on(
    bg = "auto",
    fg = "auto",
    accent = "auto",
    font = "Helvetica",
    sequential = sequential_gradient(),
    qualitative = okabe_ito(),
    inherit = FALSE
  )
  
  

  # loader for plots and tables
  w <- Waiter$new(id = c("dataSinglePlate", "ddctAbsGraph"), html = spin_loader(), color = transparent(.5))

  # Button UpdateCheck
  observeEvent(input$updateCheck, {
    runningVersion <- "0.1.3"
    url <- ("https://www.uibk.ac.at/microbiology/services/qrat/latest_version.txt")
    latestVersion <- readLines(url)
    if (runningVersion == latestVersion) {
      show_alert(
        title = "Alright,",
        text = "you're running the latest version!",
        type = "success"
      )
    } else {
      show_alert(
        title = "Update!",
        text = tags$span("A new (and better) version of qRAT is available.", tags$br(), tags$br(), tags$a(href = "https://www.uibk.ac.at/microbiology/services/qrat/", "Download here!")),
        html = TRUE,
        type = "warning"
      )
    }
  })

  # Button readPublication
  observeEvent(input$readPublication, {
    show_alert(
      title = NULL,
      text = tags$span(
        tags$h3("Publication",
          style = "color: steelblue;"
        ),
        tags$br(),
        tags$a(href = "#", "Link")
      ),
      html = TRUE
    )
  })
  
  # Button housekeepingInfo
  observeEvent(input$housekeepingInfo, {
    show_alert(
      title = "Reference Genes",
      text = tags$span("The stability of all appointed reference genes needs to be validated in advance. Popular algorithms to determine the most stable reference (housekeeping) genes 
			from a set of candidate reference are geNorm, BestKeeper and NormFinder. On average 2-4 reference genes should ideally be used for final normalization in a given experiment.",
        tags$br(),
        tags$a(href = "https://doi.org/10.1186%2Fgb-2002-3-7-research0034", "More Information")
      ),
      html = TRUE
    )
  })
  
  # Button housekeepingInfo Multiple Plates
  observeEvent(input$housekeepingInfoMP, {
    show_alert(
      title = "Reference Genes",
      text = tags$span("The stability of all appointed reference genes needs to be validated in advance. Popular algorithms to determine the most stable reference (housekeeping) genes 
			from a set of candidate reference are geNorm, BestKeeper and NormFinder. On average 2-4 reference genes should ideally be used for final normalization in a given experiment.",
                       tags$br(),
                       tags$a(href = "https://doi.org/10.1186%2Fgb-2002-3-7-research0034", "More Information")
      ),
      html = TRUE
    )
  })

  ####
  # Data Processing
  ####

  # Single Plate read file

  setData <- reactive({
    req(input$dtfile)
    f <- input$dtfile
    if (is.null(f)) {
      return(NULL)
    } else {
      filex <- f$datapath
    }
    sep <- "\t"
    if (input$Sep == "comma") sep <- ","
    if (input$Sep == "semicolon") sep <- ";"


    # Scan through beginning of file, max 100 lines and look for row with "Well"
    file.header <- readLines(con = filex, n = 100)
    n.header <- grep("^Well", file.header) - 1
    if (length(n.header) == 0) {
      n.header <- 0
    }

    xx <- fread(filex, header = TRUE, stringsAsFactor = FALSE, sep = sep, skip = n.header)

    ## change all column names to lowercase and rename columns to Well, Sample, Gene, Ct
    ## extend lists with column names from different qPCR machines in future versions of qRAT

    xx <- xx %>%
      rename_with(tolower) %>%
      rename(
        Sample = any_of(c("sample", "sample name", "name")),
        Well = any_of(c("well")),
        Gene = any_of(c("gene", "target", "target name", "detector", "primer/probe")),
        Ct = any_of(c("ct", "cq"))
      )


    ## generate replicate number for dataview
    xx$TempRepNum <- paste(xx$Sample, xx$Gene)
    xx$rp.num <- ave(xx$Sample, xx$TempRepNum, FUN = seq_along)
    xx <- subset(xx, select = c(Well, Sample, Gene, Ct, rp.num)) %>%
      filter(Well != "")

    ## replace "," with "." for dataview
    CorrectCt <- xx$Ct
    CorrectCt <- gsub(",", ".", CorrectCt, ignore.case = TRUE)
    xx$Ct <- as.numeric(CorrectCt)


    ## remove (set to NA) replicates based on threshold (bad reps and max ct) and empty rows (no name, no gene)
    correctRep <- xx %>%
      group_by(Sample, Gene) %>%
      mutate(Ct = remove_replicates(Ct)) %>%
      mutate(Ct = remove_maxCt(Ct)) %>%
      filter(Sample != "") %>%
      filter(Gene != "") %>%
      ungroup()

    badRep <- correctRep %>%
      filter(is.na(Ct)) %>%
      select(Well, Sample, Gene, rp.num)

    origData <- xx

    # for spatial plate view
    plateView <- xx %>%
      separate(Well,
        into = c("text", "num"),
        sep = "(?<=[A-Za-z])(?=[0-9])"
      )

    xx <- correctRep


    write.table(xx, file = "fileSingle.csv", sep = sep)

    Genes <- unique(xx$Gene)
    refs <- grep("^(ACT|UBQ|REF|RF).*", Genes, ignore.case = TRUE)
    if (length(refs) == 0) {
      refs <- Genes[length(Genes)]
    } else {
      refs <- Genes[refs]
    }
    # Update Inputfield for Gene Input
    updatePickerInput(session, "Refs", choices = Genes, selected = refs)

    ## read to qPCR format single plate
    htset <- NULL
    ddset <- NULL
    xdata <- read.qPCRtable("fileSingle.csv", sep = sep)
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))

      # Update Inputfields (Plot Settings and Input Settings)

      updateSelectInput(session, "Mock", choices = Samples, selected = Samples[1])
      updateSelectInput(session, "Comps", choices = Samples, selected = Samples) # update input statistics

      if (length(Samples) > 4 | max(nchar(Samples)) > 6) {
        updateSliderInput(session, "side1", value = 8)
        updateSliderInput(session, "xsrt", value = 30)
      }
    }

    return(list(data = xx, data.ht = htset, data.ddct = ddset, badReplicates = badRep, originalData = origData, PV = plateView, CompSamples = Samples))
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


  # Multiple Plate read files

  multiData <- reactive({
    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"

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

    for (i in 1:numfiles) {
      file <- ff[[i, "datapath"]]
      file.header <- readLines(con = file, n = 100)
      n.header <- grep("^Well", file.header) - 1
      if (length(n.header) == 0) {
        n.header <- 0
      }
      filelist[[i]] <- fread(file, header = TRUE, stringsAsFactor = FALSE, sep = sep, skip = n.header)
    }


    # multiplePlates <- rbindlist(lapply(f, fread, skip=n.header, sep=sep), use.names = TRUE, fill = TRUE, idcol = "PlateNumber")

    multiplePlates <- rbindlist(filelist, use.names = TRUE, fill = TRUE, idcol = "PlateNumber")


    ## change all column names to lowercase and rename columns to Well, Sample, Gene, Ct
    ## extend lists with column names from different qPCR machines
    multiplePlates <- multiplePlates %>%
      rename_with(tolower) %>%
      rename(
        Sample = any_of(c("sample", "sample name", "name")),
        Well = any_of(c("well")),
        Gene = any_of(c("gene", "target", "target name", "detector", "primer/probe")),
        Ct = any_of(c("ct", "cq")),
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
    updateSelectInput(session, "PlateSelect", choices = ChosenPlate, selected = ChosenPlate[1])

    # prepare data for spatial plate view
    MultiplePlatesView <- multiplePlates %>%
      separate(Well,
        into = c("text", "num"),
        sep = "(?<=[A-Za-z])(?=[0-9])"
      )

    ## remove (set to NA) replicates based on threshold (bad reps and max ct) and empty rows (no name, no gene)
    correctRep <- multiplePlates %>%
      group_by(PlateNumber, Sample, Gene) %>%
      mutate(Ct = remove_replicatesMP(Ct)) %>%
      mutate(Ct = remove_maxCtMulti(Ct)) %>%
      filter(Sample != "") %>%
      filter(Gene != "") %>%
      ungroup()

    badRep <- correctRep %>%
      filter(is.na(Ct)) %>%
      select(PlateNumber, Well, Sample, Gene, rp.num)

    origData <- multiplePlates

    multiplePlates <- correctRep



    ## read Samples for IPC extraction
    multiplePlatesCOPY <- multiplePlates
    SamplesWithIPCs <- unique(as.character(multiplePlatesCOPY$Sample))
    updateSelectInput(session, "IPC", choices = SamplesWithIPCs, selected = SamplesWithIPCs[0])



    write.table(multiplePlates, file = "fileMulti.csv", sep = sep)


    return(list(data = multiplePlates, originalData = origData, badReplicates = badRep, MPV = MultiplePlatesView))
  })


  ## multiData 2, for separating IPCs from main Samples

  multiData2 <- reactive({
    req(input$plates)
    # req(input$IPC)

    if (input$SepM == "tab") sep <- "\t"
    if (input$SepM == "comma") sep <- ","
    if (input$SepM == "semicolon") sep <- ";"

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

      write.table(info2, file = "fileMultiNoIPCs.csv", sep = sep)
    } else {

      # don't try to remove IPCs from samples
      write.table(info2, file = "fileMultiNoIPCs.csv", sep = sep)
    }


    ## read genes
    Genes <- unique(info2$Gene)
    refs <- grep("^(ACT|UBQ|REF|RF).*", Genes, ignore.case = TRUE)
    if (length(refs) == 0) {
      refs <- Genes[length(Genes)]
    } else {
      refs <- Genes[refs]
    }
    updatePickerInput(session, "RefsM", choices = Genes, selected = refs)



    ## read to qPCR format multiple plates
    htset <- NULL
    ddset <- NULL
    xdata <- read.qPCRtableMulti("fileMultiNoIPCs.csv", sep = sep)
    if (!is.null(xdata)) {
      htset <- xdata$data.ht
      ddset <- xdata$data.ddct
      Samples <- unique(as.character(pData(htset)$Sample))

      # Update Inputfields (Sample Selection)

      updateSelectInput(session, "MockM", choices = Samples, selected = Samples[1])
      updateSelectInput(session, "CompsM", choices = Samples, selected = Samples) # update input statistics

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

    # Update Inputfields (Plot Settings and Input Settings)
    updatePickerInput(session, "SamplePicker", choices = SamplesdCt, selected = SamplesdCt)
    updatePickerInput(session, "SamplePickerDDCt", choices = SamplesddCt, selected = SamplesddCt)

    #prepare ddCt output for statistical analysis
    ddCt_stat <- expr.rel %>% select(Sample, Gene, ddCt)
    ddCt_stat <- rename(ddCt_stat, Ct = ddCt)
    
    
    write.table(ddCt_stat, file = "ddCtsingle.csv")
    
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

    # Update Inputfields (Plot Settings and Input Settings)
    updatePickerInput(session, "SamplePickerMulti", choices = SamplesdCtMulti, selected = SamplesdCtMulti)
    updatePickerInput(session, "SamplePickerDDCtMulti", choices = SamplesddCtMulti, selected = SamplesddCtMulti)

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

      extractedIPCs <- na.omit(extractedIPCs)
      write.table(extractedIPCs, file = "IPCs.csv", sep = sep)

      gm_mean <- function(x, na.rm = TRUE) {
        exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
      }

      geometricMean <- gm_mean(extractedIPCs$Ct)

      extractedIPCs <- extractedIPCs %>%
        mutate(CtDivGmean = Ct / geometricMean)

      CalibrationFactors <- extractedIPCs %>%
        group_by(PlateNumber) %>%
        dplyr::summarize(Mean = mean(CtDivGmean, na.rm = TRUE))


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
      updatePickerInput(session, "SamplePickerIPCcomparison", choices = SamplesIPCcomparison, selected = "")
      
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
            action = DT::JS(
              "function ( e, dt, node, config ) {
          dt.page.len(-1);
          dt.ajax.reload();}"
            )
          ),
          list(
            extend = "",
            text = "Show Less",
            action = DT::JS(
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
      responsive = FALSE,
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


  output$dataSinglePlate <- DT::renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- setData()

    if (is.null(info)) {
      NULL
    } else {
      setData()$originalData
    }
  })

  ## Raw Data Table Plate Plot Spatial Original Data

  output$plotCtCard <- renderPlotly({
    info <- setData()
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
    info <- setData()

    if (is.null(info)) NULL

    fig <- ggplot(info$PV, aes(x = Ct, colour = Sample), size=0.6, scale_color_viridis()) +
      geom_density() +
      geom_rug() +
      scale_color_brewer(palette = "Set1")+
      theme_pubr(base_size=4.76 * .pt)+
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

  output$dataSinglePlateBadRep <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- setData()
    if (is.null(info)) NULL else setData()$badReplicates
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
        xaxis = list(linewidth = 2, zeroline = FALSE, ticks = "outside")
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

  output$ddctAbsolute <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
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
        PlotDataPick <- "-dCq"
        yTitle <- "-\u0394Cq"
      } else {
        PlotDataError <- "dCt.sd"
        PlotDataPick <- "dCt"
        yTitle <- "\u0394Cq"
      }
    }

    df <- info$absolute %>%
      filter(Sample %in% input$SamplePicker)
    
    
    if (PlotType == "Bar Chart") {
      figNormal <- plot_ly(df[order(df$Gene), ],
                           x = ~Sample, y = ~ get(PlotDataPick), color = ~Gene, type = "bar", error_y = list(array = ~ get(PlotDataError), color = "#000000"),
                           colors = colorPick, marker = list(size = 10, line = list(color = "rgba(0, 0, 0, .8)", width = 2))) %>%
        layout(barmode = "group", bargroupgap = 0.1)
      
      if (scalePick == "normal") {
        fig <- figNormal
      } else {
        fig <- figNormal %>% layout(yaxis = list(type = "log", dtick = 1))
      }
      
    } else {
      
      
      figNormal <- ggplot(df[order(df$Gene), ], aes(Sample, get(PlotDataPick), 
                                                    ymin = get(PlotDataPick) - get(PlotDataError), 
                                                    ymax = get(PlotDataPick) + get(PlotDataError)))+
        geom_point(
          aes(fill = Gene),
          position = position_dodge(0.6), size = 3, color = "black", shape = 21, stroke = 0.25)+
        geom_errorbar(aes(fill = Gene), width=.7, color = "black", position = position_dodge(0.6))+
        scale_fill_brewer(palette = colorPick)+
        theme_pubr(base_size=4.76 * .pt)+
        list(ggplottheme)
      
      
      if (scalePick == "normal") {
        fig <- ggplotly(figNormal)
      } else {
        fig <- figNormal +
          scale_y_log10(labels = scales::comma_format(big.mark = ""))
        fig <- ggplotly(fig)
      }
      
    }
    
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = TRUE, displayModeBar = TRUE
      ) %>%
      layout(xaxis = xaxis, yaxis = yaxis, font=t, plot_bgcolor = "transparent", margin = list(t = 100, l = 100)) %>%
      layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t, legend = list(title = list(text = "Genes")))
    
    fig
  })

  ## ddCt Expression Data Table Single Plate

  output$ddctRelative <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
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

    SamplePicker <- input$SamplePickerDDCt
    df <- info$relative %>% filter(Sample %in% c(SamplePicker))

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
    


    if (PlotType == "Bar Chart") {
      figNormal <- plot_ly(df[order(df$Gene), ],
                           x = ~Sample, y = ~ get(PlotDataPick), color = ~Gene, type = "bar", error_y = list(array = ~ get(PlotDataError), color = "#000000"),
                           colors = colorPick, marker = list(size = 10, line = list(color = "rgba(0, 0, 0, .8)", width = 2))) %>%
        layout(barmode = "group", bargroupgap = 0.1)
      
      if (scalePick == "normal") {
        fig <- figNormal
      } else {
        fig <- figNormal %>% layout(yaxis = list(type = "log", dtick = 1))
      }
      
    } else {
      
      
      figNormal <- ggplot(df[order(df$Gene), ], aes(Sample, get(PlotDataPick), 
                                                    ymin = get(PlotDataPick) - get(PlotDataError), 
                                                    ymax = get(PlotDataPick) + get(PlotDataError)))+
        geom_point(
          aes(fill = Gene),
          position = position_dodge(0.6), size = 3, color = "black", shape = 21, stroke = 0.25)+
        geom_errorbar(aes(fill = Gene), width=.7, color = "black", position = position_dodge(0.6))+
        scale_fill_brewer(palette = colorPick)+
        theme_pubr(base_size=4.76 * .pt)+
        list(ggplottheme)
      
      
      if (scalePick == "normal") {
        fig <- ggplotly(figNormal)
      } else {
        fig <- figNormal +
          scale_y_log10(labels = scales::comma_format(big.mark = ""))
        fig <- ggplotly(fig)
      }
      
    }
    
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = TRUE, displayModeBar = TRUE
      ) %>%
      layout(xaxis = xaxis, yaxis = yaxis, font=t, plot_bgcolor = "transparent", margin = list(t = 100, l = 100)) %>%
      layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t, legend = list(title = list(text = "Genes")))
    
    fig
  })



  ## MULTIPLE PLATES OUTPUT


  ## Raw Data Table Multiple Plates

  output$multiplePlatesData <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- multiData()
    if (is.null(info)) NULL else info$originalData %>% rename(Plate = PlateNumber)
  })


  ## Raw Data Table Plate Plot Spatial Original Data Multiple Plates

  output$MultiplotCtCard <- renderPlotly({
    info <- multiData()
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
      theme_pubr(base_size=4.76 * .pt)+
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

  output$dataMultiplePlatesBadRep <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
    info <- multiData()
    if (is.null(info)) NULL else multiData()$badReplicates %>% rename(Plate = PlateNumber)
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
        xaxis = list(linewidth = 2, zeroline = FALSE, ticks = "outside")
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

  output$extractedIPCsTable <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
    if (input$Id027 == TRUE) {
      info <- extractIPCs()
      if (is.null(info)) NULL else select(extractIPCs()$extrIPCs, -TempRepNum) # remove TempRepNum for displaying the table because this was just an internal variable
    } else {}
  })



  ## Calibration Factors Table

  output$tableCalibrationFactors <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, {
    if (input$Id027 == TRUE) {
      data <- extractIPCs()
      if (is.null(data)) NULL else data$CalibrationFactors %>% mutate(across(where(is.numeric), round, 4))
    } else {}
  })

  ## Inter-Plate Ct values, calibrated with calibration factor

  output$interPlateCalibration <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
    if (input$Id027 == TRUE) {
      data <- interPlateCalibration()
      if (is.null(data)) NULL else data %>% mutate(across(where(is.numeric), round, 2))
    } else {}
  })


  ## Comparison of calibrated data with non-calibrated data

  output$comparisonCalibrated <- renderPlotly({
    if (input$Id027 == TRUE) {
      info <- interPlateCalibration_Plot()
      if (is.null(info)) NULL

      calibrated <- info$calibrated  %>% 
        group_by(Sample_Gene, PlateNumber) %>%
        dplyr::summarize(MeanCq = mean(Ct, na.rm=TRUE))
      calibrated_sd <- info$calibrated  %>% 
        group_by(Sample_Gene, PlateNumber) %>%
        dplyr::summarize(sdCq = sd(Ct, na.rm = TRUE))
      calibrated <- data.frame(calibrated, calibrated_sd$sdCq)
      calibrated <- rename(calibrated, c("sdCq" = "calibrated_sd.sdCq"))
      
      
      noncalibrated <- info$noncalibrated  %>% 
        group_by(Sample_Gene, PlateNumber) %>%
        dplyr::summarize(MeanCq = mean(Ct, na.rm=TRUE))
      noncalibrated_sd <- info$noncalibrated  %>% 
        group_by(Sample_Gene, PlateNumber) %>%
        dplyr::summarize(sdCq = sd(Ct, na.rm = TRUE))
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

  output$ddctAbsoluteMulti <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
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
    PlotDataPick <- input$PlotDataMulti
    PlotType <- input$PlotTypeMulti
    PlotWidth <- input$width_dCqMulti
    PlotHeight <- input$height_dCqMulti
    PlotScale <- input$scale_dCqMulti


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


    df <- info$multiAbsolute %>% filter(Sample %in% c(SamplePicker))


    if (PlotType == "Bar Chart") {
      figNormal <- plot_ly(df[order(df$Gene), ],
                           x = ~Sample, y = ~ get(PlotDataPick), color = ~Gene, type = "bar", error_y = list(array = ~ get(PlotDataError), color = "#000000"),
                           colors = colorPick, marker = list(size = 10, line = list(color = "rgba(0, 0, 0, .8)", width = 2))) %>%
        layout(barmode = "group", bargroupgap = 0.1)
      
      if (scalePick == "normal") {
        fig <- figNormal
      } else {
        fig <- figNormal %>% layout(yaxis = list(type = "log", dtick = 1))
      }
      
    } else {
      
      
      figNormal <- ggplot(df[order(df$Gene), ], aes(Sample, get(PlotDataPick), 
                                                    ymin = get(PlotDataPick) - get(PlotDataError), 
                                                    ymax = get(PlotDataPick) + get(PlotDataError)))+
        geom_point(
          aes(fill = Gene),
          position = position_dodge(0.6), size = 3, color = "black", shape = 21, stroke = 0.25)+
        geom_errorbar(aes(fill = Gene), width=.7, color = "black", position = position_dodge(0.6))+
        scale_fill_brewer(palette = colorPick)+
        theme_pubr(base_size=4.76 * .pt)+
        list(ggplottheme)
      
      
      if (scalePick == "normal") {
        fig <- ggplotly(figNormal)
      } else {
        fig <- figNormal +
          scale_y_log10(labels = scales::comma_format(big.mark = ""))
        fig <- ggplotly(fig)
      }
      
    }
    
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = TRUE, displayModeBar = TRUE
      ) %>%
      layout(xaxis = xaxis, yaxis = yaxis, font=t, plot_bgcolor = "transparent", margin = list(t = 100, l = 100)) %>%
      layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t, legend = list(title = list(text = "Genes")))
    
    fig
  })

  ## ddCq Gene Expression Table Multiple Plates

  output$ddctRelativeMulti <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
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
    PlotDataPick <- input$PlotDataDDCtMulti
    PlotType <- input$PlotTypeDDCtMulti
    PlotWidth <- input$width_ddCqMulti
    PlotHeight <- input$height_ddCqMulti
    PlotScale <- input$scale_ddCqMulti

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

    df <- info$multiRelative %>% filter(Sample %in% c(SamplePicker))

    if (PlotType == "Bar Chart") {
      figNormal <- plot_ly(df[order(df$Gene), ],
                           x = ~Sample, y = ~ get(PlotDataPick), color = ~Gene, type = "bar", error_y = list(array = ~ get(PlotDataError), color = "#000000"),
                           colors = colorPick, marker = list(size = 10, line = list(color = "rgba(0, 0, 0, .8)", width = 2))) %>%
        layout(barmode = "group", bargroupgap = 0.1)
      
      if (scalePick == "normal") {
        fig <- figNormal
      } else {
        fig <- figNormal %>% layout(yaxis = list(type = "log", dtick = 1))
      }
      
    } else {
      
      
      figNormal <- ggplot(df[order(df$Gene), ], aes(Sample, get(PlotDataPick), 
                                                    ymin = get(PlotDataPick) - get(PlotDataError), 
                                                    ymax = get(PlotDataPick) + get(PlotDataError)))+
        geom_point(
          aes(fill = Gene),
          position = position_dodge(0.6), size = 3, color = "black", shape = 21, stroke = 0.25)+
        geom_errorbar(aes(fill = Gene), width=.7, color = "black", position = position_dodge(0.6))+
        scale_fill_brewer(palette = colorPick)+
        theme_pubr(base_size=4.76 * .pt)+
        list(ggplottheme)
      
      
      if (scalePick == "normal") {
        fig <- ggplotly(figNormal)
      } else {
        fig <- figNormal +
          scale_y_log10(labels = scales::comma_format(big.mark = ""))
        fig <- ggplotly(fig)
      }
      
    }
    
    fig <- fig %>%
      config(
        displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"),
        toImageButtonOptions = list(format = formatChoice, scale = PlotScale, width=PlotWidth, height=PlotHeight), scrollZoom = TRUE, displayModeBar = TRUE
      ) %>%
      layout(xaxis = xaxis, yaxis = yaxis, font=t, plot_bgcolor = "transparent", margin = list(t = 100, l = 100)) %>%
      layout(yaxis = list(title = yTitle), xaxis = list(title = "Sample"), font=t, legend = list(title = list(text = "Genes")))
    
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
    }

    output$countTable <- renderTable(countReplicates)

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
              div(class = "float-right", icon("times", "fa-3x")),
              h4(class = "alert-heading", "Replicates")
            ),
            p(div(class = "border border-warning", strong("Following samples are below recommended number of technical replicates:"), tableOutput("countTable")))
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

    output$countTable <- renderTable(countReplicates)

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
              div(class = "float-right", icon("times", "fa-3x")),
              h4(class = "alert-heading", "Replicates")
            ),
            p(div(class = "border border-warning", strong("Following samples are below recommended number of technical replicates:"), tableOutput("countTable")))
          )
        })
      )
    }
  })



  ## Limma Outputs

  output$resultLimma <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
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

  output$resultLimmaMulti <- renderDataTable(extensions = c("Buttons", "Responsive"), options = table_options(), rownames = FALSE, filter = "top", {
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
    return(!is.null(setData()))
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
          selectizeInput(
            inputId = paste0("CompsMultiple", i),
            label = paste0("Comparison ", i),
            choices = info$CompSamples,
            multiple = TRUE,
            options = list(maxItems = 2),
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
          selectizeInput(
            inputId = paste0("CompsMultipleM", i),
            label = paste0("Comparison ", i),
            choices = info$CompSamples,
            multiple = TRUE,
            options = list(maxItems = 2),
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




  ## running version Information of qRAT for output
  output$runningVersion <- renderText({
    runningVersion <- "0.1.3"
  })
}
