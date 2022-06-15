patrick::with_parameters_test_that(
  "test prof.to.features",
  {
    testdata <- file.path("..", "testdata", .test_name)
    input_path <- file.path(testdata, filename)
    extracted_features <- readRDS(input_path)

    actual <- prof.to.features(
      extracted_features,
      sd.cut = sd_cut,
      sigma.ratio.lim = sigma_ratio_lim,
      do.plot = FALSE
    )

    expected_path <- file.path(testdata, expected_filename)
    expected <- readRDS(expected_path)
    expect_equal(actual, expected)
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0_cdf.Rds"),
      expected_filename = "mbr_test0_features.Rds",
      sd_cut = c(0.1, 100),
      sigma_ratio_lim = c(0.1, 10)
    ),
    RCX_01_shortened_v2 = list(
      filename = c("RCX_01_shortened_v2_cdf.Rds"),
      expected_filename = "RCX_01_shortened_v2_features.Rds",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100)
    )
  )
)