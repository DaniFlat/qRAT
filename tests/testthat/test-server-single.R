# Single plate workflow against the real server function (shiny::testServer).

test_that("results are complete: values without SD are kept, column n, SEM from real n", {
  quietly(testServer(server, {
    setup_single(session)
    # exclude 3 of the 4 ref G replicates of sample C (wells A1-A3) -> only 1 valid replicate, no SD
    session$setInputs(plotCtCard_selected = c("1", "2", "3"), apply_plate_selection = 1, filterMode = "manual")
    d <- resultDdct()
    abs <- d$absolute
    rel <- d$relative

    # sample C: values are shown although the reference has no SD
    expect_equal(sum(abs$Sample == "C"), 7)
    expect_true(all(is.na(abs$dCt.sd[abs$Sample == "C"])))
    expect_true(any(rel$Sample == "T" & rel$Gene == "GOI 2"))

    valid <- setData()$data %>% dplyr::filter(!is.na(Ct)) %>% dplyr::count(Sample, Gene)
    chk <- merge(abs[, c("Sample", "Gene", "n")], valid, by = c("Sample", "Gene"), suffixes = c("", ".expected"))
    expect_equal(nrow(chk), nrow(abs))
    expect_equal(chk$n, chk$n.expected)

    session$setInputs(SamplePicker = unique(as.character(abs$Sample)), GenePicker = unique(as.character(abs$Gene)),
                      PlotData = "RQ", error_type = "sem")
    pd <- prepared_dCq_data()
    expect_equal(pd$plot_df$sem_temp, pd$plot_df$expr.sd / sqrt(pd$plot_df$n))
  }))
})

test_that("automatic filter removes the replicate furthest from the mean and keeps at least two", {
  quietly(testServer(server, {
    setup_single(session)
    # C / ref G: 21.08, 20.02, 19.11, 19.55 -> 21.08 removed, then 20.02 -> 19.11 and 19.55 remain
    ct <- setData()$data %>% dplyr::filter(Sample == "C", Gene == "ref G") %>% dplyr::pull(Ct)
    expect_equal(sort(ct[!is.na(ct)]), c(19.11, 19.55))
  }))
})

test_that("user selections survive a filter change (no re-sent input updates)", {
  quietly(testServer(server, {
    # record which input updates are actually sent; testServer only sees a snapshot of server
    # variables, so wrap send_if_changed in the live server environment
    srv_env <- environment(send_if_changed)
    orig_send <- srv_env$send_if_changed
    updates <- character(0)
    assign("send_if_changed", function(key, choices, update) {
      if (isTRUE(orig_send(key, choices, update))) updates <<- c(updates, key)
    }, envir = srv_env)

    setup_single(session)
    invisible(resultDdct())
    expect_true(all(c("Refs", "Mock_Comps", "SamplePickerDDCt") %in% updates))

    session$setInputs(Refs = c("ref G", "GOI 7"), SamplePickerDDCt = "T")
    invisible(resultDdct())
    updates <- character(0)
    session$setInputs(repDev = 0.5)
    invisible(resultDdct())
    expect_false(any(c("Refs", "Mock_Comps", "SamplePickerDDCt") %in% updates))
  }))
})

test_that("p-value adjustment uses the selected method", {
  quietly(testServer(server, {
    setup_single(session)
    expect_equal(limma_adjust_method("BH"), "BH")
    expect_equal(limma_adjust_method("Holm"), "holm")
    expect_equal(limma_adjust_method("Bonferroni"), "bonferroni")
    expect_equal(limma_adjust_method(""), "BH")

    adj <- sapply(c("BH", "Holm", "Bonferroni"), function(m) {
      session$setInputs(adjustMethod = m)
      calLimma()$adj.p.value
    })
    expect_false(identical(adj[, "BH"], adj[, "Holm"]))
    expect_true(all(adj[, "BH"] <= adj[, "Holm"] + 1e-12 & adj[, "Holm"] <= adj[, "Bonferroni"] + 1e-12, na.rm = TRUE))
  }))
})

test_that("limma runs for the plots only when significance labels are shown", {
  quietly(testServer(server, {
    setup_single(session, show_significanceDDCt = FALSE)
    session$setInputs(SamplePickerDDCt = c("C", "T"), GenePickerDDCt = unique(as.character(resultDdct()$relative$Gene)))
    expect_null(prepared_ddCq_data()$stats)
    session$setInputs(show_significanceDDCt = TRUE, ShowErrorBarDDCt = TRUE)
    expect_false(is.null(prepared_ddCq_data()$stats))
    expect_match(paste(output$ddctRelGraph2, collapse = ""), "\\*\\*\\*")
  }))
})

test_that("all single-plate outputs render", {
  quietly(testServer(server, {
    setup_single(session)
    session$setInputs(SamplePicker = "T", GenePicker = unique(as.character(resultDdct()$absolute$Gene)),
                      SamplePickerDDCt = c("C", "T"), GenePickerDDCt = unique(as.character(resultDdct()$relative$Gene)))
    for (o in c("dataSinglePlate", "plotCtCard", "plotCtDistrib", "dataSinglePlateBadRep", "SinglePlateFilterPlot",
                "SinglePlateBoxplot", "SinglePlateMIQEcheck", "countTable", "ddctAbsolute", "ddctRelative",
                "dcqPlot", "ddctRelGraph2", "resultLimma")) {
      expect_no_error(output[[o]], message = o)
    }
    expect_equal(nrow(lowReplicatesSP()), 0)   # the SD filter always keeps at least two replicates
    session$setInputs(plotCtCard_selected = c("1", "2", "3"), apply_plate_selection = 1, filterMode = "manual")
    expect_equal(nrow(lowReplicatesSP()), 1)   # C / ref G with one replicate left
    expect_no_error(output$SinglePlateMIQEcheck)
  }))
})

test_that("plate view: clicked wells become manual exclusions", {
  quietly(testServer(server, {
    setup_single(session)
    session$setInputs(plotCtCard_selected = c("3", "9"), apply_plate_selection = 1)
    expect_equal(as.integer(excluded_rows()), c(3L, 9L))

    session$setInputs(filterMode = "manual")
    ct <- setData()$data$Ct
    expect_true(all(is.na(ct[c(3, 9)])))
    expect_equal(sum(is.na(ct)), 2)
    expect_match(paste(output$plotCtCard, collapse = ""), '"selected":\\["3","9"\\]')
  }))
})

test_that("removing a sample group removes exactly that group", {
  quietly(testServer(server, {
    setup_single(session)
    session$setInputs(open_grouping_modal = 1)
    for (i in 1:3) session$setInputs(add_manual_group_btn = i)
    session$setInputs(remove_grp_btn_2 = 1)
    expect_equal(active_group_ids(), c("1", "3"))
  }))
})
