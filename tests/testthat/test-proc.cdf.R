patrick::with_parameters_test_that(
  "test proc.cdf",
  {
    if(ci_skip == TRUE) skip_on_ci()
    
    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "input", filename)
    
    sut <- proc.cdf(
      input_path,
      min_presence = min_presence,
      min_elution_length = min_elution_length,
      mz_tol = mz_tol,
      intensity_weighted = intensity_weighted,
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
      min_presence = 0.5,
      min_elution_length = 12,
      intensity_weighted = FALSE,
      cache = FALSE,
      ci_skip = FALSE
    ),
    RCX_06_shortened = list(
      filename = c("RCX_06_shortened.mzML"),
      mz_tol = 1e-06,
      min_presence = 0.7,
      min_elution_length = 4,
      intensity_weighted = TRUE,
      cache = FALSE,
      ci_skip = FALSE
    ),
    RCX_07_shortened = list(
      filename = c("RCX_07_shortened.mzML"),
      mz_tol = 1e-06,
      min_presence = 0.7,
      min_elution_length = 4,
      intensity_weighted = TRUE,
      cache = FALSE,
      ci_skip = TRUE
    ),
    RCX_08_shortened = list(
      filename = c("RCX_08_shortened.mzML"),
      mz_tol = 1e-06,
      min_presence = 0.7,
      min_elution_length = 4,
      intensity_weighted = TRUE,
      cache = FALSE,
      ci_skip = TRUE
    )
  )
)
