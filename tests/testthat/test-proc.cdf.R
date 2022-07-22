patrick::with_parameters_test_that(
  "test proc.cdf",
  {
    if(ci_skip == TRUE) skip_on_ci()
    
    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "input", filename)
    
    actual <- proc.cdf(
      input_path,
      min.pres = min_pres,
      min.run = min_run,
      tol = tol,
      intensity.weighted = intensity_weighted,
      cache = FALSE
    )

    expected_path <- file.path(testdata, "filtered", expected_filename)
    expected <- readRDS(expected_path)

    # exclude last column from comparison as there lies the stochastic nature
    expect_equal(actual[, 1:3], expected[, 1:3])
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0.mzml"),
      expected_filename = "mbr_test0_cdf.Rds",
      tol = 1e-05,
      min_pres = 0.5,
      min_run = 12,
      intensity_weighted = FALSE,
      ci_skip = FALSE
    ),
    RCX_01_shortened_v2 = list(
      filename = c("RCX_06_shortened.mzML"),
      expected_filename = "RCX_06_shortened_cdf.Rds",
      tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE,
      ci_skip = FALSE
    ),
    RCX_09_shortened_v2 = list(
      filename = c("RCX_07_shortened.mzML"),
      expected_filename = "RCX_07_shortened_cdf.Rds",
      tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE,
      ci_skip = TRUE
    ),
    RCX_16_shortened_v2 = list(
      filename = c("RCX_08_shortened.mzML"),
      expected_filename = "RCX_08_shortened_cdf.Rds",
      tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE,
      ci_skip = TRUE
    )
  )
)
