test_that("basic two-step hybrid test", {
    test_names <- c("mbr_test0.mzml",
                    "mbr_test1.mzml",
                    "mbr_test2.mzml",
                    "mbr_test0_copy.mzml")
    test_path <- paste0("../testdata/", test_names)
    metadata <- "../testdata/two_step_hybrid_info.csv"

    tempdir <- tempdir()
    temp_path <- paste0(tempdir, "/", test_names)
    file.copy(test_path, temp_path)

    expected_final_ftrs <- readRDS("../testdata/final_ftrs.Rda")
    known_table <- arrow::read_parquet('../testdata/known_table.parquet')

    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
    if (nzchar(chk) && chk == "TRUE") {
        num_workers <- 2L
    } else {
        num_workers <- parallel::detectCores()
    }

    result <- two.step.hybrid(
        filenames = test_names,
        metadata = metadata,
        work_dir = tempdir,
        known.table = known_table,
        cluster = num_workers)
    final_features <- result$final_features

    expected_final_features <- as.data.frame(expected_final_ftrs)
    colnames(expected_final_features)[2:4] <- c("rt", "mz_min", "mz_max")
    colnames(expected_final_features) <- stringr::str_remove_all(colnames(expected_final_features), ".mzml")
    keys <- c("mz", "rt", "mz_min", "mz_max")
    final_features <- as_tibble(arrange_at(final_features, keys))
    expected_final_features <- as_tibble(arrange_at(expected_final_features, keys))

    comparison <- dataCompareR::rCompare(
        final_features,
        expected_final_features,
        keys = keys
    )

    dataCompareR::saveReport(
        comparison,
        reportName = "final_features_comparison",
        reportLocation = "..",
        showInViewer = FALSE,
        missmatchCount = 10000
    )

    expect_equal(final_features, expected_final_features)
    expect_equal(final_features$mz, expected_final_features$mz, tolerance=0.001)

})