patrick::with_parameters_test_that(
  "test compute_clusters",
  {
    testdata <- file.path("..", "testdata")

    filenames <- sapply(files, function(x) {
      file.path(testdata, input, paste0(x, ".parquet"))
    })

    extracted <- lapply(filenames, function(x) {
      tibble::as_tibble(arrow::read_parquet(x)) |> dplyr::rename(rt = pos)
    })

    actual <- compute_clusters(
      feature_tables = extracted,
      mz_tol_relative = NA,
      rt_tol_relative = NA,
      mz_max_diff = mz_max_diff,
      mz_tol_absolute = mz_tol_absolute
    )

    expected <- lapply(files, function(x) {
      filepath <- file.path(testdata, "clusters", paste0(x, "_", input, "_clusters.parquet"))
      tibble::as_tibble(arrow::read_parquet(filepath))
    })


    for(i in seq_along(filenames)) {
      expect_equal(actual$feature_tables[[i]], expected[[i]])
    }

  },
  patrick::cases(
    RCX_shortened_extracted = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      input = "extracted",
      mz_max_diff = 10 * 1e-05,
      mz_tol_absolute = 0.01
    ),
    RCX_shortened_adjusted = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      input = "adjusted",
      mz_max_diff = 10 * 1e-05,
      mz_tol_absolute = 0.01
    )
  )
)
