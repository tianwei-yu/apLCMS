patrick::with_parameters_test_that(
  "test proc.cdf test",
  {
    testdata <- file.path("..", "testdata", .test_name)
    input_path <- file.path(testdata, filename)
    actual <- proc.cdf(
      input_path,
      min.pres = min_pres,
      min.run = min_run,
      tol = tol,
      intensity.weighted = intensity_weighted,
      cache = FALSE
    )

    expected_path <- file.path(testdata, expected_filename)
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
      intensity_weighted = FALSE
    ),
    RCX_01_shortened_v2 = list(
      filename = c("RCX_01_shortened_v2.mzML"),
      expected_filename = "RCX_01_shortened_v2_cdf.Rds",
      tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE
    ),
    RCX_09_shortened_v2 = list(
      filename = c("RCX_09_shortened_v2.mzML"),
      expected_filename = "RCX_09_shortened_v2_cdf.Rds",
      tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE
    ),
    RCX_16_shortened_v2 = list(
      filename = c("RCX_16_shortened_v2.mzML"),
      expected_filename = "RCX_16_shortened_v2_cdf.Rds",
      tol = 1e-06,
      min_pres = 0.7,
      min_run = 4,
      intensity_weighted = TRUE
    )
  )
)
