patrick::with_parameters_test_that(
  "test cont.index.R",
  {
    if(ci_skip == TRUE) skip_on_ci()

    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "input", filename)

    input_data <- as.matrix(arrow::read_parquet(input_path) ) 
    actual <- cont.index(input_data, min.pres, min.run)

    actual <- tibble::tibble(
      mz = actual$new.rec[, 1],
      rt = actual$new.rec[, 2],
      intensity = actual$new.rec[, 3],
      group_number = actual$new.rec[, 4]
      )

    expected_path <- file.path(testdata, "filtered", "cont_index", paste0(.test_name, "_cont_index.parquet"))
    expected <- arrow::read_parquet(expected_path) 
  
    expect_equal(actual, expected)
    
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0.parquet"),
      min.pres = 0.5,
      min.run = 12,
      ci_skip = FALSE
    )
  )
)
