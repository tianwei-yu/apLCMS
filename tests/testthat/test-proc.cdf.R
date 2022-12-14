patrick::with_parameters_test_that(
  "test proc.cdf",
  {
    if(ci_skip == TRUE) skip_on_ci()

    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "input", filename)

    sut <- proc.cdf(
      input_path,
      min_pres = min_pres,
      min_run = min_run,
      mz_tol = mz_tol,
      baseline_correct = 0.0,
      baseline_correct_noise_percentile = 0.05,
      intensity_weighted = intensity_weighted,
      do.plot = FALSE,
      cache = cache
    )

    expected_path <- file.path(testdata, "filtered", paste0(.test_name, ".parquet"))

    # exclude last column from comparison as there lies the stochastic nature
    expected <- arrow::read_parquet(expected_path) |> dplyr::select(-group_number)
    actual <- sut |> dplyr::select(-group_number)

    expect_equal(actual, expected)
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0.mzml"),
      mz_tol = 1e-05,
      min_pres = 0.5,
      min_run = 12,
      intensity_weighted = FALSE,
      cache = FALSE,
      ci_skip = FALSE
    ),
    RCX_06_shortened = list(
      filename = c("RCX_06_shortened.mzML"),
      mz_tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE,
      cache = FALSE,
      ci_skip = FALSE
    ),
    RCX_07_shortened = list(
      filename = c("RCX_07_shortened.mzML"),
      mz_tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE,
      cache = FALSE,
      ci_skip = TRUE
    ),
    RCX_08_shortened = list(
      filename = c("RCX_08_shortened.mzML"),
      mz_tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE,
      cache = FALSE,
      ci_skip = TRUE
    )
  )
)
