patrick::with_parameters_test_that("basic hybrid test", {
  if(ci_skip == TRUE) skip_on_ci()

  testdata <- file.path("..", "testdata")

  test_files <- sapply(files, function(x) {
    file.path(testdata, "input", x)
  })

  known_table <- arrow::read_parquet(
    file.path(testdata, "hybrid", "known_table.parquet")
  )

  actual <- hybrid(
    test_files,
    known_table,
    align_mz_tol = NA,
    align_rt_tol = NA,
    cluster = get_num_workers())

  expected <- arrow::read_parquet(
    file.path(testdata, "hybrid", paste0(.test_name, "_recovered_feature_sample_table.parquet"))
  )

  expect_equal(actual$recovered_feature_sample_table, expected)
}, patrick::cases(
  mbr = list(
    files = c("mbr_test0.mzml", "mbr_test1.mzml", "mbr_test2.mzml"),
    ci_skip = TRUE
  ),
  RCX_shortened = list(
    files = c("RCX_06_shortened.mzML", "RCX_07_shortened.mzML", "RCX_08_shortened.mzML"),
    ci_skip = FALSE
  )
))
