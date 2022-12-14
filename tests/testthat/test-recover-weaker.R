patrick::with_parameters_test_that(
  "recover weaker signals test",
  {
    store_reports <- FALSE
    testdata <- file.path("..", "testdata")

    ms_files <- sapply(files, function(x) {
      file.path(testdata, "input", paste0(x, ".mzML"))
    })

    extracted <- lapply(files, function(x) {
      xx <- file.path(testdata, "extracted", paste0(x, ".parquet"))
      arrow::read_parquet(xx)
    })

    adjusted <- lapply(files, function(x) {
      xx <- file.path(testdata, "adjusted", paste0(x, ".parquet"))
      arrow::read_parquet(xx)
    })

    aligned <- load_aligned_features(
      file.path(testdata, "aligned", "metadata_table.parquet"),
      file.path(testdata, "aligned", "intensity_table.parquet"),
      file.path(testdata, "aligned", "rt_table.parquet"),
      file.path(testdata, "aligned", "tolerances.parquet")
    )

    recovered <- lapply(seq_along(ms_files), function(i) {
      recover.weaker(
        filename = ms_files[[i]],
        sample_name = files[i],
        extracted_features = extracted[[i]],
        adjusted_features = adjusted[[i]],
        metadata_table = aligned$metadata,
        rt_table = aligned$rt,
        intensity_table = aligned$intensity,
        mz_tol = mz_tol,
        mz_tol_relative = aligned$mz_tol_relative,
        rt_tol_relative = aligned$rt_tol_relative,
        recover_mz_range = recover_mz_range,
        recover_rt_range = recover_rt_range,
        use_observed_range = use_observed_range,
        bandwidth = bandwidth,
        min_bandwidth = min_bandwidth,
        max_bandwidth = max_bandwidth,
        recover_min_count = recover_min_count,
        intensity_weighted = intensity_weighted
      )
    })

    # create and load final files
    keys <- c("mz", "rt", "sd1", "sd2", "area")

    extracted_recovered_actual <- lapply(recovered, function(x) x$extracted_features |> dplyr::arrange_at(keys))
    corrected_recovered_actual <- lapply(recovered, function(x) x$adjusted_features |> dplyr::arrange_at(keys))


    extracted_recovered_expected <- lapply(files, function(x) {
      xx <- file.path(testdata, "recovered", "recovered-extracted", paste0(x, ".parquet"))
      arrow::read_parquet(xx) |> dplyr::arrange_at(keys)
    })


    corrected_recovered_expected <- lapply(files, function(x) {
      xx <- file.path(testdata, "recovered", "recovered-corrected", paste0(x, ".parquet"))
      arrow::read_parquet(xx) |> dplyr::arrange_at(keys)
    })
    

    # compare files
    for (i in seq_along(files)) {
      # extracted recovered
      actual_extracted_i <- extracted_recovered_actual[[i]]
      expected_extracted_i <- extracted_recovered_expected[[i]]

      expect_equal(actual_extracted_i, expected_extracted_i)

      # corrected recovered
      actual_corrected_i <- corrected_recovered_actual[[i]]
      expected_corrected_i <- corrected_recovered_expected[[i]]

      expect_equal(actual_corrected_i, expected_corrected_i)

      if (store_reports) {
        report_extracted <- dataCompareR::rCompare(
          actual_extracted_i,
          expected_extracted_i,
          keys = keys,
          roundDigits = 4,
          mismatches = 100000
        )
        dataCompareR::saveReport(
          report_extracted,
          reportName = paste0(files[[i]], "_extracted"),
          showInViewer = FALSE,
          HTMLReport = FALSE,
          mismatchCount = 10000
        )
        report_corrected <- dataCompareR::rCompare(
          actual_corrected_i,
          expected_corrected_i,
          keys = keys,
          roundDigits = 4,
          mismatches = 100000
        )
        dataCompareR::saveReport(
          report_corrected,
          reportName = paste0(files[[i]], "_adjusted"),
          showInViewer = FALSE,
          HTMLReport = FALSE,
          mismatchCount = 10000
        )
      }
    }
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      mz_tol = 1e-05,
      recover_mz_range = NA,
      recover_rt_range = NA,
      use_observed_range = TRUE,
      min_bandwidth = NA,
      max_bandwidth = NA,
      recover_min_count = 3,
      bandwidth = 0.5,
      intensity_weighted = FALSE
    )
  )
)