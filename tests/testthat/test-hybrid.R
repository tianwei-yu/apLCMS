test_that("basic hybrid test", {
  test_files <- c('../testdata/mbr_test0.mzml',
                  '../testdata/mbr_test1.mzml',
                  '../testdata/mbr_test2.mzml')
  
  expected <- arrow::read_parquet('../testdata/hybrid_recovered_feature_sample_table.parquet')
  known_table <- arrow::read_parquet('../testdata/known_table.parquet')
  
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
  }

  result <- hybrid(test_files, known_table, cluster = num_workers)

  expect_equal(result$recovered_feature_sample_table, expected)
})
