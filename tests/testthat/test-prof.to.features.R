patrick::with_parameters_test_that(
  "test prof.to.features",
  {
    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "filtered", filename)
    extracted_features <- readRDS(input_path)

    actual <- prof.to.features(
      extracted_features,
      sd.cut = sd_cut,
      sigma.ratio.lim = sigma_ratio_lim,
      shape.model = shape_model,
      do.plot = FALSE
    )

    expected_path <- file.path(testdata, "features", expected_filename)
    expected <- readRDS(expected_path)
    
    expect_equal(actual, expected)
  },
  patrick::cases(
    # mbr_test0 = list(
    #   filename = c("mbr_test0_cdf.Rds"),
    #   expected_filename = "mbr_test0_features.Rds",
    #   sd_cut = c(0.1, 100),
    #   sigma_ratio_lim = c(0.1, 10),
    #   shape_model = "bi-Gaussian"
    # ),
    mbr_test0_gaussian = list(
      filename = c("mbr_test0_cdf.Rds"),
      expected_filename = "mbr_test0_gaussian_features.Rds",
      sd_cut = c(0.1, 100),
      sigma_ratio_lim = c(0.1, 10),
      shape_model = "bi-Gaussian"
    ),
    RCX_01_shortened_gaussian = list(
      filename = c("RCX_06_shortened_cdf.Rds"),
      expected_filename = "RCX_06_shortened_gaussian_features.Rds",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "Gaussian"
    ),
    RCX_01_shortened_v2 = list(
      filename = c("RCX_06_shortened_cdf.Rds"),
      expected_filename = "RCX_06_shortened_features.Rds",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "bi-Gaussian"
    ),
    RCX_09_shortened_v2 = list(
      filename = c("RCX_07_shortened_cdf.Rds"),
      expected_filename = "RCX_07_shortened_features.Rds",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "bi-Gaussian"
    ),
    RCX_16_shortened_v2 = list(
      filename = c("RCX_08_shortened_cdf.Rds"),
      expected_filename = "RCX_08_shortened_features.Rds",
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      shape_model = "bi-Gaussian"
    )
  )
)
