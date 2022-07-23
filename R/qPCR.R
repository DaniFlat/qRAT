####
# qPCR - Relative Expression Analysis Tool
# version: 0.1.5
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
# read qPCR data with HTqPCR and calculate with ddCt packages
####


## read qPCRTable Single Plate

read.qPCRtable <- function(fname, na.value = 40, ...) {
  if (missing(fname)) stop("No input fname specified")
  dt.raw <- fread(fname, header = TRUE, ...)
  for (cc in colnames(dt.raw)) {
    dcc <- dt.raw[, ..cc]
    if (is.character(dcc)) dcc[is.na(dcc)] <- "NA"
    dt.raw[, cc] <- dcc
  }


  ## check replicate column
  cnames <- colnames(dt.raw)
  tbs <- table(dt.raw$rp.num)
  if (min(tbs) != max(tbs)) showNotification("Number of replicates should be identical accross all samples!", type = "error", duration = 5)



  ## check samples and genes
  if (!any(grepl("^Sample$", cnames))) {
    dt.raw$Sample <- apply(dt.raw[, grep("^ds\\.", cnames)], 1, paste, collapse = "@")
  }
  dt.raw$lb.Sample <- paste(dt.raw$Sample, dt.raw$rp.num, sep = "@")
  tbs <- table(dt.raw$lb.Sample)
  if (min(tbs) != max(tbs)) showNotification("All samples should have identical list of genes", type = "error", duration = 5)
  tbs <- table(dt.raw$Gene)
  if (min(tbs) != max(tbs)) showNotification("Each gene should be tested in all samples", type = "error", duration = 5)

  # handle Cq values
  xct <- dt.raw$Ct

  # replace "," with "."
  xct <- gsub(",", ".", xct, ignore.case = TRUE)


  xct <- gsub("(Undetermined)|(No Ct)", "60", xct, ignore.case = TRUE)


  xct <- as.numeric(xct)
  xct[xct > na.value] <- NA
  dt.raw$Ct <- xct

  Genes <- unique(dt.raw$Gene)
  dt.raw$Gene <- factor(dt.raw$Gene, levels = Genes)
  dx <- dcast(dt.raw, Gene ~ lb.Sample, value.var = "Ct")
  X <- as.matrix(dx[, -1])
  nspots <- nrow(X)
  nSamples <- ncol(X)

  Samples <- sapply(colnames(X), function(x) {
    paste(rev(rev(strsplit(x, split = "@", fixed = TRUE)[[1]])[-1]), collapse = "@")
  })
  reps <- sapply(colnames(X), function(x) {
    rev(strsplit(x, split = "@", fixed = TRUE)[[1]])[1]
  })

  ## qPCRset object
  pdata <- data.frame(Sample = Samples, rep = reps, row.names = 1:nSamples)
  Sample.info <- new("AnnotatedDataFrame", data = pdata)
  featPos <- paste0("A", 1:nspots)
  featType <- rep("Target", nspots)
  featType[grep("act|ref|ubq", Genes, ignore.case = TRUE)] <- "Reference"
  featType <- factor(featType)

  df <- data.frame(featureNames = Genes, featureType = featType, featurePos = featPos)
  metaData <- data.frame(labelDescription = c(
    "Name of the qPCR feature (Gene)",
    "Type of feature", "Position on assay"
  ))
  featData <- AnnotatedDataFrame(data = df, varMetadata = metaData)

  colnames(X) <- NULL
  X.flags <- data.frame(matrix("Passed", ncol = nSamples, nrow = nspots))
  X.cat <- data.frame(matrix("OK", ncol = nSamples, nrow = nspots),
    stringsAsFactors = FALSE
  )
  htset <- new("qPCRset",
    exprs = X, phenoData = Sample.info, featureData = featData,
    featureCategory = X.cat, flag = X.flags, CtHistory = data.frame()
  )

  ## InputFrame object (ddCt)
  coreData <- data.frame(matrix(nrow = nrow(dt.raw), ncol = 4))
  colnames(coreData) <- c("Sample", "Detector", "Ct", "Platename")

  ddset <- new("InputFrame")
  ddset@coreData <- coreData
  ddset@coreData$Sample <- dt.raw$Sample
  ddset@coreData$Ct <- dt.raw$Ct
  ddset@coreData$Detector <- dt.raw$Gene
  ddset@coreData$Platename <- rep(1, nrow(dt.raw))
  ddset@files <- fname

  list(data.ht = htset, data.ddct = ddset)
}



## read qPCRTable Multiple Plates

read.qPCRtableMulti <- function(fname, na.value = 40, ...) {
  if (missing(fname)) stop("No input fname specified")
  dt.raw <- fread(fname, header = TRUE, ...)
  for (cc in colnames(dt.raw)) {
    dcc <- dt.raw[, ..cc]
    if (is.character(dcc)) dcc[is.na(dcc)] <- "NA"
    dt.raw[, cc] <- dcc
  }


  ## add replicate number
  dt.raw$TempRepNum <- paste(dt.raw$Sample, dt.raw$Gene)
  dt.raw$rp.num <- ave(dt.raw$Sample, dt.raw$TempRepNum, FUN = seq_along)

  ## check column names
  cnames <- colnames(dt.raw)
  tbs <- table(dt.raw$rp.num)


  ## check samples and genes
  if (!any(grepl("^Sample$", cnames))) {
    dt.raw$Sample <- apply(dt.raw[, grep("^ds\\.", cnames)], 1, paste, collapse = "@")
  }
  dt.raw$lb.Sample <- paste(dt.raw$Sample, dt.raw$rp.num, sep = "@")
  tbs <- table(dt.raw$lb.Sample)
  tbs <- table(dt.raw$Gene)

  # handle Cq values
  xct <- dt.raw$Ct
  xct <- gsub("(Undetermined)|(No Ct)", "60", xct, ignore.case = TRUE)

  # replace "," with "."
  xct <- gsub(",", ".", xct, ignore.case = TRUE)
  xct <- as.numeric(xct)
  xct[xct > na.value] <- NA
  dt.raw$Ct <- xct

  Genes <- unique(dt.raw$Gene)
  dt.raw$Gene <- factor(dt.raw$Gene, levels = Genes)
  dx <- dcast(dt.raw, Gene ~ lb.Sample, value.var = "Ct")
  X <- as.matrix(dx[, -1])
  nspots <- nrow(X)
  nSamples <- ncol(X)

  Samples <- sapply(colnames(X), function(x) {
    paste(rev(rev(strsplit(x, split = "@", fixed = TRUE)[[1]])[-1]), collapse = "@")
  })
  reps <- sapply(colnames(X), function(x) {
    rev(strsplit(x, split = "@", fixed = TRUE)[[1]])[1]
  })

  ## qPCRset object
  pdata <- data.frame(Sample = Samples, rep = reps, row.names = 1:nSamples)
  Sample.info <- new("AnnotatedDataFrame", data = pdata)
  featPos <- paste0("A", 1:nspots)
  featType <- rep("Target", nspots)
  featType[grep("act|ref|ubq", Genes, ignore.case = TRUE)] <- "Reference"
  featType <- factor(featType)

  df <- data.frame(featureNames = Genes, featureType = featType, featurePos = featPos)
  metaData <- data.frame(labelDescription = c(
    "Name of the qPCR feature (Gene)",
    "Type of feature", "Position on assay"
  ))
  featData <- AnnotatedDataFrame(data = df, varMetadata = metaData)

  colnames(X) <- NULL
  X.flags <- data.frame(matrix("Passed", ncol = nSamples, nrow = nspots))
  X.cat <- data.frame(matrix("OK", ncol = nSamples, nrow = nspots),
    stringsAsFactors = FALSE
  )
  htset <- new("qPCRset",
    exprs = X, phenoData = Sample.info, featureData = featData,
    featureCategory = X.cat, flag = X.flags, CtHistory = data.frame()
  )

  ## InputFrame object (ddCt)
  coreData <- data.frame(matrix(nrow = nrow(dt.raw), ncol = 4))
  colnames(coreData) <- c("Sample", "Detector", "Ct", "Platename")

  ddset <- new("InputFrame")
  ddset@coreData <- coreData
  ddset@coreData$Sample <- dt.raw$Sample
  ddset@coreData$Ct <- dt.raw$Ct
  ddset@coreData$Detector <- dt.raw$Gene
  ddset@coreData$Platename <- rep(1, nrow(dt.raw))
  ddset@files <- fname

  list(data.ht = htset, data.ddct = ddset)
}



#Calculate and return major results of HTqPCR

calHTqPCR <- function(dt, ref.Gene, comp.strings = NULL, comp.type, adjustMethod, norm.method = "deltaCt", verbose = FALSE) {
  results <- list(result = NULL, contrast = NULL, d.norm = NULL, g.norm = NULL)
  if (is.null(dt)) {
    return(results)
  }
  res <- NULL
  vss <- NULL
  if (!ref.Gene %in% fData(dt)$featureNames) {
    return(results)
  }
  d.norm <- normalizeCtData(dt, deltaCt.genes = ref.Gene, norm = norm.method, verbose = verbose)
  results$d.norm <- d.norm

  g.norm <- normalizeCtData(dt, norm = "geometric.mean", verbose = FALSE)
  results$g.norm <- g.norm


  if (is.null(comp.strings)) {
    return(results)
  }
  Samples <- pData(dt)$Sample
  treatments <- gsub("[^a-z0-9]", ".", Samples, ignore.case = TRUE)
  design <- model.matrix(~ 0 + treatments)
  colnames(design) <- gsub("treatments", "", colnames(design))

  vxx <- gsub("(^\\s+|\\s+$)", "", comp.strings)
  if (vxx == "") {
    return(results)
  }
  if (comp.type == 2) {
    if (!grepl(":", vxx)) {
      return(results)
    }
    sps <- strsplit(vxx, "[ ,;;,]+")[[1]]
    vxx <- sapply(sps, function(x) {
      aa <- strsplit(x, ":", fixed = TRUE)[[1]]
      aa <- gsub("[^a-z0-9_]", ".", aa, ignore.case = TRUE)
    })
    if (all(vxx %in% treatments)) {
      vss <- apply(vxx, 2, paste, collapse = " - ")
    } else {
      return(results)
    }
  } else {
    Samples <- unique(Samples)
    if (vxx %in% Samples) {
      sps <- Samples[Samples != vxx]
      trs <- gsub("[^a-z0-9_]", ".", sps, ignore.case = TRUE)
      vxx <- gsub("[^a-z0-9_]", ".", vxx, ignore.case = TRUE)
      vss <- paste(trs, vxx, sep = " - ")
    } else {
      return(results)
    }
  }
  names(vss) <- sps

  res <- NULL
  if (!is.null(vss)) {
    contrasts <- makeContrasts(contrasts = vss, levels = design)
    res <- limmaCtData(d.norm,
      design = design,
      contrasts = contrasts, ndups = 1, spacing = 1, adjust.method = adjustMethod
    )
  }

  results$result <- res
  results$contrast <- vss
  results
}


