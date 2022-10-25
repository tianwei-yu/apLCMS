patrick::with_parameters_test_that(
  "test cont.index.R",
  {
    if(ci_skip == TRUE) skip_on_ci()

    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "input", filename)
    input_data <- arrow::read_parquet(input_path) 
    tab_prt <- as.matrix(input_data) 
    actual <- cont.index(tab_prt, min.pres, min.run)

    tab_prt <- tibble::tibble(
      mz = actual$new.rec[, 1],
      rt = actual$new.rec[, 2],
      intensity = actual$new.rec[, 3],
      group_number = actual$new.rec[, 4]
      )

    expected_path <- file.path(testdata, "input", "output_cont.parquet")
    expected <- arrow::read_parquet(expected_path) 
  
    expect_equal(tab_prt, expected)
  },
  patrick::cases(
    input_cont = list(
      filename = c("input_cont.parquet"),
      min.pres = 0.5,
      min.run = 12,
      ci_skip = FALSE
    )
  )
)
