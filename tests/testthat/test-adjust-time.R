patrick::with_parameters_test_that(
  "get_template",
  {
    testdata <- file.path("..", "testdata")

    filenames <- lapply(files, function(x) {
      file.path(testdata, "clusters", paste0(x, "_extracted_clusters.parquet"))
    })

    extracted <- lapply(filenames, arrow::read_parquet)
    template_features <- compute_template(extracted)

    expected <- file.path(testdata, "template", "RCX_shortened.parquet")
    expected <- arrow::read_parquet(expected) |> dplyr::rename(sample_id = label)

    expect_equal(template_features, expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened")
    )
  )
)

patrick::with_parameters_test_that(
  "correct time test",
  {
    testdata <- file.path("..", "testdata")

    template_features <- file.path(testdata, "template", "RCX_shortened.parquet")
    template_features <- arrow::read_parquet(template_features) |> dplyr::rename(sample_id = label)

    extracted <- file.path(testdata, "clusters", paste0(.test_name, "_extracted_clusters.parquet"))
    extracted <- arrow::read_parquet(extracted)

    corrected <- correct_time(
      this.feature = extracted,
      template_features = template_features,
      mz_tol_relative = mz_tol_relative,
      rt_tol_relative = rt_tol_relative
    )

    expected <- tibble::as_tibble(arrow::read_parquet(file.path(testdata, "adjusted", paste0(.test_name, ".parquet"))))
    expected <- expected |> dplyr::rename( rt = pos, sample_id = V6, cluster = V7)
    expect_equal(corrected, expected)
  },
  patrick::cases(
    RCX_06_shortened = list(
      mz_tol_relative = 6.856763e-06,
      rt_tol_relative = 3.618581
    ),
    RCX_07_shortened = list(
      mz_tol_relative = 6.856763e-06,
      rt_tol_relative = 3.618581
    ),
    RCX_08_shortened = list(
      mz_tol_relative = 6.856763e-06,
      rt_tol_relative = 3.618581
    )
  )
)
