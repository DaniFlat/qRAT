#' read qPCR data with HTqPCR and calculate with ddCt packages

# data.table is attached via library() but not imported in NAMESPACE. Without this flag,
# `dt[, (cc) := ...]` below fails inside the package namespace (pkgload::load_all / installed package)
# with "[ was called on a data.table in an environment that is not data.table-aware".
.datatable.aware <- TRUE


## Input helper: accepts a file path or an in-memory data.frame.
## In-memory input avoids writing shared temp files into the working directory
## (concurrent sessions would otherwise overwrite each other's data).
read_qPCR_input <- function(fname, ...) {
  if (!is.data.frame(fname)) {
    return(fread(fname, header = TRUE, ...))
  }
  dt.raw <- as.data.table(fname)
  # mimic the former fwrite/fread round trip: factors became character, missing strings became ""
  for (cc in colnames(dt.raw)) {
    dcc <- dt.raw[[cc]]
    if (is.factor(dcc)) dcc <- as.character(dcc)
    if (is.character(dcc)) {
      dcc[is.na(dcc)] <- ""
      set(dt.raw, j = cc, value = dcc)
    }
  }
  dt.raw
}


## read qPCRTable Single Plate

read.qPCRtable <- function(fname, na.value = 40, ...) {
  if (missing(fname)) stop("No input fname specified")
  dt.raw <- read_qPCR_input(fname, ...)
  for (cc in colnames(dt.raw)) {
    dcc <- dt.raw[[cc]] 
    
    if (is.character(dcc)) {
      dcc[is.na(dcc)] <- "NA"
      dt.raw[, (cc) := dcc] 
    }
  }
  
  
  ## check replicate column
  cnames <- colnames(dt.raw)
  tbs <- table(dt.raw$rp.num)
  if (min(tbs) != max(tbs)) showNotification("Number of replicates should be identical accross all samples!", type = "error", duration = 5)
  
  
  
  ## check samples and genes
  if (!any(grepl("^Sample$", cnames))) {
    dt.raw$Sample <- apply(dt.raw[, grep("^ds\\.", cnames), with = FALSE], 1, paste, collapse = "@")
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
  ddset@files <- if (is.character(fname)) fname else "in-memory"
  
  list(data.ht = htset, data.ddct = ddset)
}



## read qPCRTable Multiple Plates

read.qPCRtableMulti <- function(fname, na.value = 40, ...) {
  if (missing(fname)) stop("No input fname specified")
  dt.raw <- read_qPCR_input(fname, ...)
  for (cc in colnames(dt.raw)) {
    dcc <- dt.raw[[cc]] 
    
    if (is.character(dcc)) {
      dcc[is.na(dcc)] <- "NA"
      dt.raw[, (cc) := dcc] 
    }
  }
  
  
  ## add replicate number
  dt.raw$TempRepNum <- paste(dt.raw$Sample, dt.raw$Gene)
  dt.raw$rp.num <- ave(dt.raw$Sample, dt.raw$TempRepNum, FUN = seq_along)
  
  ## check column names
  cnames <- colnames(dt.raw)
  tbs <- table(dt.raw$rp.num)
  
  
  ## check samples and genes
  if (!any(grepl("^Sample$", cnames))) {
    dt.raw$Sample <- apply(dt.raw[, grep("^ds\\.", cnames), with = FALSE], 1, paste, collapse = "@")
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
  ddset@files <- if (is.character(fname)) fname else "in-memory"
  
  list(data.ht = htset, data.ddct = ddset)
}



## deltaCt normalisation per biological sample (used instead of HTqPCR::normalizeCtData(norm = "deltaCt")).
## HTqPCR subtracts the reference Cq of the same replicate column (replicate 1 - replicate 1, ...): when the
## replicate filter had removed that reference replicate, the target replicate became NA and limma could end
## up without any data for a sample/gene. Here every target replicate is normalised against the mean of all
## valid reference Cq values of its sample - the same basis the ddCq calculation (ddCt package) uses.
normalize_deltaCt_by_sample <- function(dt, ref.Gene) {
  # explicit Biobase:: - server.R attaches rlang after HTqPCR, which masks exprs()
  X <- Biobase::exprs(dt)
  ref_rows <- Biobase::featureNames(dt) %in% ref.Gene
  samples <- as.character(Biobase::pData(dt)$Sample)
  ref_mean <- vapply(samples, function(s) mean(X[ref_rows, samples == s], na.rm = TRUE), numeric(1))
  ref_mean[is.nan(ref_mean)] <- NA
  Biobase::exprs(dt) <- sweep(X, 2, ref_mean, "-")
  dt
}

## Automatic replicate filter (user guide, "Automatic Filtering"): while the standard deviation of the
## technical replicates exceeds max_sd, the value furthest from their mean is removed. At least min_n
## replicates are kept - with two discordant values it cannot be decided which one is the outlier.
## Returns x with the removed values set to NA.
filter_replicates_sd <- function(x, max_sd, min_n = 2) {
  y <- x
  repeat {
    ok <- which(!is.na(y))
    if (length(ok) <= min_n || sd(y[ok]) <= max_sd) break
    y[ok[which.max(abs(y[ok] - mean(y[ok])))]] <- NA
  }
  y
}

## Biological grouping (user guide, "Biological Replicates and Grouping"): the technical replicates of each
## biological sample are averaged first, so every biological replicate enters its group with one value and
## mean, SD/SEM and the statistics are based on the biological replicates.
## d needs the columns Sample (group name), BioSample (original sample name), Gene, Ct and Well.
aggregate_biological <- function(d) {
  d %>%
    dplyr::filter(!is.na(Ct)) %>%
    dplyr::group_by(Sample, BioSample, Gene) %>%
    dplyr::summarise(Well = dplyr::first(Well), Ct = mean(Ct), .groups = "drop") %>%
    dplyr::group_by(Sample, Gene) %>%
    dplyr::mutate(rp.num = as.character(dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

## limma statistics for the qPCR data: the core of HTqPCR::limmaCtData (lmFit -> contrasts.fit -> eBayes ->
## topTable per contrast) without its per-group mean columns. Those used rowMeans() and failed as soon as a
## sample had a single column (e.g. one biological replicate after grouping); qRAT never used them.
## Returns one data.frame per contrast with genes, t.test, p.value, adj.p.value, ddCt and FC.
limma_qpcr <- function(q, design, contrasts, adjust.method) {
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(Biobase::exprs(q), design = design), contrasts = contrasts))
  out <- lapply(colnames(contrasts), function(coef) {
    res <- limma::topTable(fit, coef = coef, number = nrow(fit), sort.by = "none", adjust.method = adjust.method)
    data.frame(genes = rownames(res), t.test = res$t, p.value = res$P.Value, adj.p.value = res$adj.P.Val,
               ddCt = res$logFC, FC = 2^(-res$logFC), row.names = rownames(res), stringsAsFactors = FALSE)
  })
  names(out) <- colnames(contrasts)
  out
}

## Significance symbols for adjusted p-values, shared by the statistics table and the plots:
## ns (> 0.05), * (<= 0.05), ** (<= 0.01), *** (<= 0.001), **** (<= 0.0001); NA -> ""
signif_symbol <- function(p) {
  out <- as.character(cut(p, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                          labels = c("****", "***", "**", "*", "ns"), right = TRUE))
  out[is.na(out)] <- ""
  out
}


## Sample name -> valid R name for the limma design/contrasts ("24h" -> "X24h", "WT T1" -> "WT.T1").
## Previously names starting with a digit made makeContrasts() fail (and with it the plot with significance labels).
limma_name <- function(x) make.names(gsub("[^A-Za-z0-9_]", ".", x))


#Calculate and return major results of HTqPCR

calHTqPCR <- function(dt, ref.Gene, comp.strings = NULL, comp.type, adjustMethod, norm.method = "deltaCt", verbose = FALSE) {
  results <- list(result = NULL, contrast = NULL, d.norm = NULL, g.norm = NULL)
  if (is.null(dt) || length(ref.Gene) == 0) {
    return(results)
  }
  res <- NULL
  vss <- NULL

  # (norm.method is kept for compatibility; the app always normalises with deltaCt)
  d.norm <- normalize_deltaCt_by_sample(dt, ref.Gene)
  results$d.norm <- d.norm
  
  # g.norm (geometric mean normalisation) is not used anywhere in the app, so it is no longer computed
  
  
  if (is.null(comp.strings)) {
    return(results)
  }
  Samples <- pData(dt)$Sample
  # same naming rule for design and comparison strings
  treatments <- limma_name(Samples)
  # a comparison needs at least two samples (one sample: model.matrix/makeContrasts would fail)
  if (length(unique(treatments)) < 2) {
    return(results)
  }
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
      aa <- limma_name(aa)
    })
    if (all(vxx %in% treatments)) {
      vss <- apply(vxx, 2, paste, collapse = " - ")
      # readable labels with the original sample names, e.g. "WT T1 vs WT C" instead of "WT.T1:WT.C"
      orig <- unique(as.character(Samples))
      clean <- limma_name(orig)
      sps <- apply(vxx, 2, function(pair) paste(orig[match(pair, clean)], collapse = " vs "))
    } else {
      return(results)
    }
  } else {
    Samples <- unique(Samples)
    if (vxx %in% Samples) {
      sps <- Samples[Samples != vxx]
      trs <- limma_name(sps)
      vxx <- limma_name(vxx)
      vss <- paste(trs, vxx, sep = " - ")
    } else {
      return(results)
    }
  }
  names(vss) <- sps
  
  res <- NULL
  if (!is.null(vss)) {
    contrasts <- makeContrasts(contrasts = vss, levels = design)
    res <- limma_qpcr(d.norm, design = design, contrasts = contrasts, adjust.method = adjustMethod)
  }
  
  results$result <- res
  results$contrast <- vss
  results
}


