patrick::with_parameters_test_that(
  "recover weaker signals test",
  {
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
        sample_name = as.character(i),
        extracted_features = extracted[[i]],
        adjusted_features = adjusted[[i]],
        metadata_table = aligned$metadata,
        rt_table = aligned$rt,
        intensity_table = aligned$intensity,
        orig.tol = mz_tol,
        align.mz.tol = aligned$mz_tolerance,
        align.rt.tol = aligned$rt_tolerance,
        recover_mz_range = recover_mz_range,
        recover_rt_range = recover_rt_range,
        use.observed.range = use_observed_range,
        bandwidth = bandwidth,
        min.bw = min_bandwidth,
        max.bw = max_bandwidth,
        recover.min.count = recover_min_count,
        intensity.weighted = intensity.weighted
      )
    })

    # feature_table <- dplyr::select(aligned$metadata, c(mz, rt, mzmin, mzmax))
    # rt_crosstab <- cbind(feature_table, sapply(recovered, function(x) x$this.times))
    # int_crosstab <- cbind(feature_table, sapply(recovered, function(x) x$this.ftrs))

    # feature_names <- rownames(feature_table)
    # sample_names <- get_sample_name(ms_files)

    # browser()
    # recovered_actual <- list(
    #   extracted_features = lapply(recovered, function(x) x$this.f1),
    #   corrected_features = lapply(recovered, function(x) x$this.f2),
    #   rt = as_feature_crosstab(feature_names, sample_names, rt_crosstab),
    #   intensity = as_feature_crosstab(feature_names, sample_names, int_crosstab)
    # )

    # aligned_feature_sample_table_actual <- create_feature_sample_table(aligned)
    # recovered_feature_sample_table_actual <- create_feature_sample_table(recovered_actual)

    # aligned_feature_sample_table_expected <- arrow::read_parquet(file.path(testdata, "recovered", "aligned_feature_sample_table.parquet"))
    # recovered_feature_sample_table_expected <- arrow::read_parquet(file.path(testdata, "recovered", "recovered_feature_sample_table.parquet"))

    # expect_equal(aligned_feature_sample_table_actual, aligned_feature_sample_table_expected)
    # expect_equal(recovered_feature_sample_table_actual, recovered_feature_sample_table_expected)

    # create and load final files

    extracted_recovered_actual <- lapply(recovered, function(x) x$extracted_features)
    corrected_recovered_actual <- lapply(recovered, function(x) x$adjusted_features)

    filenames <- lapply(files, function(x) {
      file.path(testdata, "recovered", "recovered-extracted", paste0(x, ".parquet"))
    })

    extracted_recovered_expected <- lapply(filenames, arrow::read_parquet)

    filenames <- lapply(files, function(x) {
      file.path(testdata, "recovered", "recovered-corrected", paste0(x, ".parquet"))
    })

    corrected_recovered_expected <- lapply(filenames, function(x) {
       arrow::read_parquet(x)
    })
    # preprocess dataframes
    keys <- c("mz", "rt", "sd1", "sd2", "area")

    extracted_recovered_actual <- lapply(extracted_recovered_actual, function(x) {
      x |> dplyr::arrange_at(keys)
    })
    corrected_recovered_actual <- lapply(corrected_recovered_actual, function(x) {
      x |> dplyr::arrange_at(keys)
    })

    extracted_recovered_expected <- lapply(extracted_recovered_expected, function(x) {
      x |> dplyr::arrange_at(keys)
    })
    corrected_recovered_expected <- lapply(corrected_recovered_expected, function(x) {
      x |> dplyr::arrange_at(keys)
    })

    # compare files
    for (i in seq_along(files)) {
      # extracted recovered
      actual_extracted_i <- extracted_recovered_actual[[i]]
      expected_extracted_i <- extracted_recovered_expected[[i]]

      expect_equal(actual_extracted_i, expected_extracted_i)

      # report <- dataCompareR::rCompare(actual_extracted_i, expected_extracted_i, keys = keys, roundDigits = 4, mismatches = 100000)
      # dataCompareR::saveReport(report, reportName = paste0(files[[i]],"_extracted"), showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

      # corrected recovered
      actual_corrected_i <- corrected_recovered_actual[[i]]
      expected_corrected_i <- corrected_recovered_expected[[i]]

      expect_equal(actual_corrected_i, expected_corrected_i)

      # report <- dataCompareR::rCompare(actual_corrected_i, expected_corrected_i, keys = keys, roundDigits = 4, mismatches = 100000)
      # dataCompareR::saveReport(report, reportName = paste0(files[[i]],"_adjusted"), showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)
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
      intensity.weighted = FALSE
    )
  )
)