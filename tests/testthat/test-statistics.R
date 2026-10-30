# Statistics helpers (qPCR.R) and the limma normalisation.

test_that("significance symbols use one scheme for table and plots", {
  p <- c(0, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.03, 0.05, 0.051, 1, NA)
  expect_equal(signif_symbol(p),
               c("****", "****", "****", "***", "***", "**", "**", "*", "*", "ns", "ns", ""))
})

test_that("deltaCt uses the mean reference Cq of each sample", {
  quietly(testServer(server, {
    setup_multi(session)
    ht <- multiData2()$data.ht
    norm <- normalize_deltaCt_by_sample(ht, "ref G")
    X <- Biobase::exprs(ht); N <- Biobase::exprs(norm)
    s <- Biobase::pData(ht)$Sample; g <- Biobase::featureNames(ht)
    cols <- which(s == "WT T1")
    ref_mean <- mean(X[g == "ref G", cols], na.rm = TRUE)
    expect_equal(unname(N[g == "GOI 1", cols]), unname(X[g == "GOI 1", cols] - ref_mean))
  }))
})

test_that("deltaCt has results where the matching reference replicate was filtered out", {
  # per-replicate pairing (HTqPCR) gave NA wherever the reference replicate in the same position was missing
  quietly(testServer(server, {
    setup_multi(session)
    ht <- multiData2()$data.ht
    s <- Biobase::pData(ht)$Sample; g <- Biobase::featureNames(ht)
    col <- which(s == "WT T1")[1]
    X <- Biobase::exprs(ht); X[g == "ref G", col] <- NA; Biobase::exprs(ht) <- X
    N <- Biobase::exprs(normalize_deltaCt_by_sample(ht, "ref G"))
    expect_false(is.na(N[g == "GOI 1", col]))
  }))
})

test_that("replicate filter: furthest outlier first, iterative, at least two kept", {
  expect_equal(filter_replicates_sd(c(20, 25), 0.3), c(20, 25))                  # two values: undecidable, kept
  expect_equal(filter_replicates_sd(c(20, 20.1, 22), 0.3), c(20, 20.1, NA))
  expect_equal(filter_replicates_sd(c(21.08, 20.02, 19.11, 19.55), 0.3), c(NA, NA, 19.11, 19.55))
  expect_equal(filter_replicates_sd(c(20, 20.1, 20.2), 0.3), c(20, 20.1, 20.2))  # within threshold
  expect_equal(filter_replicates_sd(c(NA, 20, 23, 20.1), 0.3), c(NA, 20, NA, 20.1))
  expect_equal(filter_replicates_sd(c(NA, NA, NA), 0.3), c(NA, NA, NA))
})

test_that("biological aggregation: one mean value per biological sample", {
  d <- data.frame(Well = paste0("A", 1:7), Sample = "grp", BioSample = rep(c("s1", "s2", "s3"), c(3, 2, 2)),
                  Gene = "g", Ct = c(20, 21, NA, 22, 24, 25, 25))
  a <- aggregate_biological(d)
  expect_equal(nrow(a), 3)
  expect_equal(a$Ct, c(20.5, 23, 25))
  expect_equal(a$rp.num, c("1", "2", "3"))
})

test_that("limma works when a sample has a single value", {
  quietly(testServer(server, {
    setup_multi(session)
    ht <- multiData2()$data.ht
    tr <- limma_name(Biobase::pData(ht)$Sample)
    keep <- tr != "WT.T1" | !duplicated(tr)   # WT T1 with one column only
    q <- ht[, keep]; tr <- tr[keep]
    design <- model.matrix(~ 0 + tr); colnames(design) <- sub("^tr", "", colnames(design))
    res <- limma_qpcr(normalize_deltaCt_by_sample(q, "ref G"), design,
                      limma::makeContrasts(contrasts = "WT.T1 - WT.C", levels = design), "BH")
    expect_false(all(is.na(res[["WT.T1 - WT.C"]]$p.value)))
  }))
})

test_that("geNorm ranks the most stable pair together on rank 1", {
  quietly(testServer(server, {
    set.seed(2)
    base <- rnorm(6, 20)
    ct <- cbind(A = base + rnorm(6, 0, 0.05), B = base + 1 + rnorm(6, 0, 0.05), C = base + rnorm(6, 0, 0.5), D = rnorm(6, 25, 2))
    ct <- rbind(ct, c(NA, 20, 20, 20))   # incomplete sample is dropped
    res <- run_genorm_logic(ct)
    expect_setequal(res$Gene[res$Rank == 1], c("A", "B"))
    expect_equal(res$Gene[res$Rank == 3], "D")
    expect_true(all(diff(res$M) >= -1e-12))
    expect_error(run_genorm_logic(ct[, 1, drop = FALSE]))
  }))
})

test_that("paired comparisons are labelled with the original sample names", {
  quietly(testServer(server, {
    setup_multi(session, compTypeM = "Multiple paired Comparisons")
    session$setInputs(add_btnM = 1, CompsMultipleM1 = c("M T1", "M C"))
    expect_equal(unique(calLimmaMulti()$Label), "M T1 vs M C")
  }))
})
