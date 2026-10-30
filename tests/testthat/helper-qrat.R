# Shared helpers for the qRAT tests.
# Test data: tests/testdata/single_plate.csv (one plate, samples C/T, 8 genes) and
# tests/testdata/plate_1.csv + plate_2.csv (two plates incl. an IPC sample "IPC").

testdata <- function(...) testthat::test_path("..", "testdata", ...)

# data.frame in the shape of a fileInput() value
upload <- function(...) {
  f <- c(...)
  data.frame(name = basename(f), size = file.size(f), type = "text/csv", datapath = f)
}

# set default inputs; arguments in ... override the defaults
set_inputs <- function(session, defaults, ...) {
  do.call(session$setInputs, utils::modifyList(defaults, list(...)))
}

# start a single-plate session with the test plate loaded
setup_single <- function(session, ...) {
  set_inputs(session, list(Sep = "auto", filterMode = "auto", maxCt = c(5, 35), repDev = 0.3, NTC_Input = NULL,
                           groupSamples = FALSE, compType = "Single Comparison", limma_input = "dCq", adjustMethod = "BH"), ...)
  session$setInputs(dtfile = upload(testdata("single_plate.csv")))
  session$setInputs(Refs = "ref G", Mock = "C", Comps = "C")
}

# start a multi-plate session with both test plates loaded
setup_multi <- function(session, ...) {
  set_inputs(session, list(SepM = "auto", filterModeMP = "auto", maxCtMulti = c(5, 35), repDevMulti = 0.3, NTC_Input_MP = NULL,
                           groupSamplesMulti = FALSE, compTypeM = "Single Comparison", limma_inputMulti = "dCq",
                           adjustMethodMulti = "BH", Id027 = FALSE, IPC = "IPC"), ...)
  session$setInputs(plates = upload(testdata("plate_1.csv"), testdata("plate_2.csv")))
  session$setInputs(RefsM = "ref G", MockM = "WT C", CompsM = "WT C")
}

# the server's analysis code emits many harmless warnings (limma, ddCt, ggplot2)
quietly <- function(expr) suppressMessages(suppressWarnings(expr))
