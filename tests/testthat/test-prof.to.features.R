patrick::with_parameters_test_that(
  "test prof.to.features",
  {
    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "filtered", filename)
    extracted_features <- arrow::read_parquet(input_path)

    actual <- prof.to.features(
      profile = extracted_features,
      sd.cut = sd_cut,
      sigma.ratio.lim = sigma_ratio_lim,
      shape.model = shape_model,
      do.plot = do.plot
    )

    expected_path <- file.path(testdata, "features", expected_filename)
    expected <- arrow::read_parquet(expected_path)
    expect_equal(actual, expected)
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0.parquet"),
      expected_filename = "mbr_test0_features.parquet",
      sd_cut = c(0.1, 100),
      sigma_ratio_lim = c(0.1, 10),
      shape_model = "bi-Gaussian",
      do.plot = FALSE
    ),
    RCX_06_shortened_gaussian = list(
      filename = c("RCX_06_shortened.parquet"),
      expected_filename = "RCX_06_shortened_gaussian_features.parquet",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "Gaussian",
      do.plot = FALSE
    ),
    RCX_06_shortened_v2 = list(
      filename = c("RCX_06_shortened.parquet"),
      expected_filename = "RCX_06_shortened_features.parquet",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "bi-Gaussian",
      do.plot = FALSE
    ),
    RCX_07_shortened_v2 = list(
      filename = c("RCX_07_shortened.parquet"),
      expected_filename = "RCX_07_shortened_features.parquet",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "bi-Gaussian",
      do.plot = FALSE
    ),
    RCX_8_shortened_v2 = list(
      filename = c("RCX_08_shortened.parquet"),
      expected_filename = "RCX_08_shortened_features.parquet",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "bi-Gaussian",
      do.plot = FALSE
    )
  )
)
