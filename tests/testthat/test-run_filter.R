patrick::with_parameters_test_that(
  "test run_filter",
  {
    if(ci_skip == TRUE) skip_on_ci()

    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "input", filename)

    input_data <- as.matrix(arrow::read_parquet(input_path) ) 
    actual <- run_filter(input_data, min_pres, min_run)

    actual <- tibble::tibble(
      mz = actual$new_rec[, 1],
      rt = actual$new_rec[, 2],
      intensity = actual$new_rec[, 3],
      group_number = actual$new_rec[, 4]
      )

    expected_path <- file.path(testdata, "filtered", "run_filter", paste0(.test_name, "_run_filter.parquet"))
    expected <- arrow::read_parquet(expected_path) 
  
    expect_equal(actual, expected)
    
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0.parquet"),
      min_pres = 0.5,
      min_run = 12,
      ci_skip = FALSE
    )
  )
)
