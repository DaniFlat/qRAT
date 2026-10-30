# Faulty or unusual uploads must not end the session (an observer error closes a Shiny session)
# and should still give results where possible.

write_lines_tmp <- function(lines, ext = "csv") {
  f <- tempfile(fileext = paste0(".", ext))
  writeLines(lines, f)
  f
}
base_lines <- function() readLines(testdata("single_plate.csv"))

upload_single <- function(session, f) {
  session$setInputs(Sep = "auto", filterMode = "auto", maxCt = c(5, 35), repDev = 0.3, NTC_Input = NULL,
                    groupSamples = FALSE, compType = "Single Comparison", limma_input = "dCq", adjustMethod = "BH")
  session$setInputs(dtfile = data.frame(name = basename(f), size = file.size(f), type = "", datapath = f))
}

test_that("upload without Cq column: no crash, no analysis", {
  f <- write_lines_tmp(sub(",[^,]*$", "", base_lines()))   # drop last column (Cq)
  quietly(testServer(server, {
    upload_single(session, f)
    expect_false(session$isClosed())
    expect_null(readSinglePlateData())
  }))
})

test_that("wrong file type (.xlsx): no crash, no analysis", {
  f <- tempfile(fileext = ".xlsx")
  writeBin(as.raw(c(0x50, 0x4b, 3, 4, rep(0, 200))), f)
  quietly(testServer(server, {
    upload_single(session, f)
    expect_false(session$isClosed())
    expect_null(readSinglePlateData())
  }))
})

test_that("file with 'Sample' and 'Sample Name' columns is read (first match wins)", {
  l <- base_lines()
  l <- paste0(l, ",", c("Sample", sub("^[^,]*,([^,]*),.*$", "\\1", l[-1])))
  f <- write_lines_tmp(l)
  quietly(testServer(server, {
    upload_single(session, f)
    expect_false(session$isClosed())
    expect_equal(sort(unique(readSinglePlateData()$data$Sample)), c("C", "T"))
  }))
})

test_that("instrument header lines before the table are skipped", {
  f <- write_lines_tmp(c("* Block Type = 96-Well Block", "* Experiment = test", "", base_lines()))
  quietly(testServer(server, {
    upload_single(session, f)
    expect_equal(nrow(readSinglePlateData()$data), 64)
  }))
})

test_that("sample names starting with a digit work in limma and keep their names", {
  l <- base_lines()
  l[-1] <- sub("^([^,]*),C,", "\\1,0h,", sub("^([^,]*),T,", "\\1,24h,", l[-1]))
  f <- write_lines_tmp(l)
  quietly(testServer(server, {
    upload_single(session, f)
    session$setInputs(Refs = "ref G", Mock = "0h", Comps = "0h", show_significanceDDCt = TRUE, ShowErrorBarDDCt = TRUE)
    res <- calLimma()
    expect_false(is.null(res))
    expect_equal(unique(res$Label), "24h")
    session$setInputs(SamplePickerDDCt = c("0h", "24h"), GenePickerDDCt = unique(as.character(resultDdct()$relative$Gene)))
    expect_match(paste(output$ddctRelGraph2, collapse = ""), "\\*\\*\\*")
  }))
})

test_that("only one sample: no statistics, but the plot still renders", {
  l <- base_lines()
  f <- write_lines_tmp(c(l[1], grep(",T,", l, value = TRUE)))
  quietly(testServer(server, {
    upload_single(session, f)
    session$setInputs(Refs = "ref G", Mock = "T", Comps = "T", show_significanceDDCt = TRUE)
    expect_null(calLimma())
    session$setInputs(SamplePickerDDCt = "T", GenePickerDDCt = unique(as.character(resultDdct()$relative$Gene)))
    expect_no_error(output$ddctRelGraph2)
    expect_error(output$resultLimma, "No statistics available")
  }))
})

test_that("multi plate: a plate without Cq column is excluded instead of analysed as empty", {
  p2 <- readLines(testdata("plate_2.csv"))
  f2 <- write_lines_tmp(sub(",[^,]*$", "", p2))
  quietly(testServer(server, {
    session$setInputs(SepM = "auto", filterModeMP = "auto", maxCtMulti = c(5, 35), repDevMulti = 0.3, NTC_Input_MP = NULL,
                      groupSamplesMulti = FALSE, compTypeM = "Single Comparison", limma_inputMulti = "dCq",
                      adjustMethodMulti = "BH", Id027 = FALSE, IPC = "IPC")
    session$setInputs(plates = upload(testdata("plate_1.csv"), f2))
    expect_false(session$isClosed())
    d <- readMPData()$data
    expect_equal(unique(d$PlateNumber), 1L)
    expect_equal(sum(is.na(d$Ct)), 0)
  }))
})
