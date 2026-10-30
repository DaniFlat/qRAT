# Multi plate workflow against the real server function (shiny::testServer).

test_that("IPC calibration feeds ddCq, limma and geNorm", {
  quietly(testServer(server, {
    setup_multi(session, Id027 = TRUE, adjustMethodMulti = "Holm")
    invisible(resultDdctMulti())
    expect_equal(multiData2()$SamplesAnalysed$Ct, interPlateCalibration()$Ct)
    expect_false(is.null(calLimmaMulti()))   # uses adjustMethodMulti, single-plate input is unset

    session$setInputs(ref_candidates_multi = c("ref G", "GOI 1", "GOI 2"), run_val_btn_multi = 1)
    expect_equal(nrow(validationResultsMulti()), 3)
  }))
})

test_that("calibration shifts ddCq by the calibration factor", {
  ddcq <- new.env()
  for (ipc in c(FALSE, TRUE)) quietly(testServer(server, {
    setup_multi(session, Id027 = ipc)
    ddcq[[as.character(ipc)]] <- resultDdctMulti()$multiRelative
  }))
  raw <- ddcq[["FALSE"]]; cal <- ddcq[["TRUE"]]
  cal <- cal[match(paste(raw$Sample, raw$Gene), paste(cal$Sample, cal$Gene)), ]
  shift <- cal$ddCt - raw$ddCt
  # test data: plate offsets are +/- 0.015, calibrator WT C lies on both plates
  expect_true(all(abs(shift) <= 0.0151, na.rm = TRUE))
  expect_true(any(abs(shift) > 0.01, na.rm = TRUE))
})

test_that("paired comparisons react to the selection in the comparison boxes", {
  quietly(testServer(server, {
    setup_multi(session, compTypeM = "Multiple paired Comparisons")
    session$setInputs(add_btnM = 1, CompsMultipleM1 = c("WT T1", "WT C"))
    expect_equal(unique(calLimmaMulti()$Label), "WT T1 vs WT C")
    session$setInputs(CompsMultipleM1 = c("WT T2", "WT C"))
    expect_equal(unique(calLimmaMulti()$Label), "WT T2 vs WT C")
    session$setInputs(CompsMultipleM1 = "WT T3")   # incomplete pair: no result, no error
    expect_null(calLimmaMulti())
  }))
})

test_that("grouping: n and SD are based on the biological replicates", {
  quietly(testServer(server, {
    setup_multi(session, groupSamplesMulti = TRUE)
    confirmedGroupingMulti(list(`1` = list(name = "WT early", samples = c("WT T1", "WT T2", "WT T3"))))
    a <- resultDdctMulti()$multiAbsolute
    expect_true(all(a$n[a$Sample == "WT early"] == 3))
    # each biological sample enters with the mean of its technical replicates
    tech <- multiData2()$SamplesCalibrated %>%
      dplyr::filter(!is.na(Ct), BioSample %in% c("WT T1", "WT T2", "WT T3"), Gene == "GOI 1") %>%
      dplyr::group_by(BioSample) %>% dplyr::summarise(Ct = mean(Ct))
    bio <- multiData2()$SamplesAnalysed %>% dplyr::filter(Sample == "WT early", Gene == "GOI 1")
    expect_equal(sort(bio$Ct), sort(tech$Ct))
    expect_false(is.null(calLimmaMulti()))
  }))
})

test_that("column n counts valid replicates across plates", {
  quietly(testServer(server, {
    setup_multi(session)
    a <- resultDdctMulti()$multiAbsolute
    valid <- multiData2()$SamplesAnalysed %>% dplyr::filter(!is.na(Ct)) %>% dplyr::count(Sample, Gene)
    chk <- merge(a[, c("Sample", "Gene", "n")], valid, by = c("Sample", "Gene"), suffixes = c("", ".expected"))
    expect_equal(nrow(chk), nrow(a))
    expect_equal(chk$n, chk$n.expected)
  }))
})

test_that("all multi-plate outputs render and plots are ggiraph widgets with clean legends", {
  quietly(testServer(server, {
    setup_multi(session, Id027 = TRUE, PlateSelect = "1", SamplePickerIPCcomparison = c("WT C_ref G", "M C_ref G"))
    a <- resultDdctMulti()
    session$setInputs(SamplePickerMulti = unique(as.character(a$multiAbsolute$Sample)), GenePickerMulti = unique(as.character(a$multiAbsolute$Gene)),
                      SamplePickerDDCtMulti = unique(as.character(a$multiRelative$Sample)), GenePickerDDCtMulti = unique(as.character(a$multiRelative$Gene)),
                      ref_candidates_multi = c("ref G", "GOI 1", "GOI 2"), run_val_btn_multi = 1)
    for (o in c("multiplePlatesData", "MultiplotCtCard", "MultiplotCtDistrib", "dataMultiplePlatesBadRep", "MultiplePlatesFilterPlot",
                "MultiplePlatesBoxplot", "extractedIPCsTable", "tableCalibrationFactors", "interPlateCalibration", "plotIPCcomparison",
                "ddctAbsoluteMulti", "ddctRelativeMulti", "ddctMultiAbsGraph", "ddctRelGraphMulti", "resultLimmaMulti",
                "countTableMulti", "ref_ranking_table_multi", "ref_stability_plot_multi")) {
      expect_no_error(output[[o]], message = o)
    }
    ipc <- jsonlite::fromJSON(output$plotIPCcomparison)
    expect_match(ipc$x$html, "<svg")
    expect_false(grepl("(Raw,1)", ipc$x$html, fixed = TRUE))
  }))
})

test_that("plot export uses the multi-plate settings", {
  quietly(testServer(server, {
    setup_multi(session, ShowErrorBarMulti = TRUE, PlotDataMulti = "RQ")
    a <- resultDdctMulti()$multiAbsolute
    session$setInputs(SamplePickerMulti = unique(as.character(a$Sample)), GenePickerMulti = unique(as.character(a$Gene)),
                      plotTitleMulti = "MULTI-TITLE", plotTitle = "SINGLE-TITLE",
                      plot_format_Multi = "png", plot_w_Multi = 6, plot_h_Multi = 4, plot_res_Multi = 1)
    f <- output$download_plot_MultiAbs
    expect_true(file.exists(f) && file.size(f) > 1000)
    d <- prepared_MultiAbs_data()
    p <- render_qPCR_plot(d$plot_df, d$stats, d$col_val, d$err_col, d$ylab, input_type = "Multi", export_mode = TRUE)
    expect_equal(p$labels$title, "MULTI-TITLE")
  }))
})

test_that("multi plate view: a selection only replaces the exclusions of the plate shown", {
  quietly(testServer(server, {
    setup_multi(session, PlateSelect = "1")
    session$setInputs(MultiplotCtCard_selected = "2", apply_plate_selection_MP = 1)
    session$setInputs(PlateSelect = "2", MultiplotCtCard_selected = c("100", "101", "5"), apply_plate_selection_MP = 2)
    expect_equal(as.integer(excluded_rows_MP()), c(2L, 100L, 101L))
  }))
})
