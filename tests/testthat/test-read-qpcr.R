# read.qPCRtable()/read.qPCRtableMulti() accept a data.frame (in memory) and give the same
# objects as the former fwrite -> fread round trip through a temp file.

as_processed <- function(f) {
  d <- as.data.frame(data.table::fread(f))
  d$rp.num <- as.character(d$rp.num)
  d
}

compare_paths <- function(fun, processed) {
  tmp <- tempfile(fileext = ".csv")
  data.table::fwrite(processed, tmp, sep = "\t")
  from_file <- quietly(fun(tmp, sep = "\t"))
  in_memory <- quietly(fun(processed))
  expect_identical(Biobase::exprs(in_memory$data.ht), Biobase::exprs(from_file$data.ht))
  expect_identical(Biobase::pData(in_memory$data.ht), Biobase::pData(from_file$data.ht))
  expect_identical(in_memory$data.ddct@coreData, from_file$data.ddct@coreData)
}

test_that("single plate: in-memory input equals file round trip", {
  compare_paths(read.qPCRtable, as_processed(testdata("legacy", "fileSingle.csv")))
})

test_that("multi plate: in-memory input equals file round trip", {
  compare_paths(read.qPCRtableMulti, as_processed(testdata("legacy", "fileMultiNoIPCs.csv")))
})
