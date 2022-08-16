get_num_workers <- function() {
  # CRAN limits the number of cores available to packages to 2
  # source https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions#50571533
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
  }
  return(num_workers)
}

test_that("basic hybrid test", {
  test_files <- c('../testdata/mbr_test0.mzml',
                  '../testdata/mbr_test1.mzml',
                  '../testdata/mbr_test2.mzml')
  
  expected <- arrow::read_parquet('../testdata/hybrid_recovered_feature_sample_table.parquet')
  known_table <- arrow::read_parquet('../testdata/known_table.parquet')

  actual <- hybrid(test_files, known_table, cluster = get_num_workers())

  actual$recovered_feature_sample_table <- actual$recovered_feature_sample_table |> dplyr::arrange_all()
  expected <- expected |> dplyr::arrange_all()

  expect_equal(actual$recovered_feature_sample_table, expected)
})
