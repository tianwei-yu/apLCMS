patrick::with_parameters_test_that(
  "feature.align test",
  {
    testdata <- file.path("..", "testdata")

    ms_files <- lapply(files, function(x) {
      file.path(testdata, "input", paste0(x, ".mzML"))
    })

    sample_names <- get_sample_name(unlist(ms_files))

    corrected_files <- lapply(files, function(x) {
      file.path(testdata, "adjusted", paste0(x, ".parquet"))
    })

    corrected_features <- lapply(corrected_files, function(x) {
      arrow::read_parquet(x)
    })

    aligned_actual <- feature.align(
      features = corrected_features,
      sample_names = files,
      min_occurrence = min_occurrence,
      mz_tol_relative = mz_tol_relative,
      rt_tol_relative = rt_tol_relative,
      mz_max_diff = 10 * mz_tol,
      mz_tol_absolute = mz_tol_absolute,
      do.plot = do.plot
    )

    aligned_expected <- load_aligned_features(
      file.path(testdata, "aligned", "metadata_table.parquet"),
      file.path(testdata, "aligned", "intensity_table.parquet"),
      file.path(testdata, "aligned", "rt_table.parquet"),
      file.path(testdata, "aligned", "tolerances.parquet")
    )

    expect_equal(aligned_actual, aligned_expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      min_occurrence = 2,
      mz_tol_relative = NA,
      rt_tol_relative = NA,
      mz_tol = 1e-05,
      mz_tol_absolute = 0.01,
      do.plot = FALSE
    )
  )
)


patrick::with_parameters_test_that(
  "compute_aligned_feature_table test",
  {
    testdata <- file.path("..", "testdata")

    ms_files <- lapply(files, function(x) {
      file.path(testdata, "clusters", paste0(x, "_adjusted_clusters.parquet"))
    })

    corrected_features <- lapply(ms_files, function(x) {
      arrow::read_parquet(x)
    })

    aligned_actual <- create_aligned_feature_table(
        dplyr::bind_rows(corrected_features),
        min_occurrence,
        files,
        rt_tol_relative,
        mz_tol_relative
    )

    aligned_expected <- list(
      metadata = arrow::read_parquet(file.path(testdata, "aligned", "metadata_table.parquet")),
      intensity = arrow::read_parquet(file.path(testdata, "aligned", "intensity_table.parquet")),
      rt = arrow::read_parquet(file.path(testdata, "aligned", "rt_table.parquet"))
    )

    expect_equal(aligned_actual, aligned_expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      min_occurrence = 2,
      mz_tol_relative = 6.85676325338646e-06,
      rt_tol_relative = 2.17918873407775
    )
  )
)
