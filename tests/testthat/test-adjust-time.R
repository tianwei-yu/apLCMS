patrick::with_parameters_test_that(
  "adjust time test",
  {
    testdata <- file.path("..", "testdata")

    filenames <- lapply(files, function(x) {
      file.path(testdata, "clusters", paste0(x, "_extracted_clusters.parquet"))
    })

    extracted <- lapply(filenames, arrow::read_parquet)

    corrected <- adjust.time(
      extracted_features = extracted,
      mz_tol_relative = mz_tol_relative,
      rt_tol_relative = rt_tol_relative,
      do.plot = do.plot
    )

    expected_filenames <- lapply(files, function(x) {
      file.path(testdata, "adjusted", paste0(x, ".parquet"))
    })

    expected <- lapply(expected_filenames, function(x) {
      tibble::as_tibble(arrow::read_parquet(x)) |> dplyr::rename(rt = pos, sample_id = V6, cluster = V7)
    })

    corrected <- lapply(corrected, tibble::as_tibble)

    expect_equal(corrected, expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      mz_tol_relative = 6.85676325338646e-06,
      rt_tol_relative = 3.61858118506494,
      do.plot = FALSE
    )
  )
)
