test_that("basic unsupervised test", {
  test_files <- c('../testdata/input/mbr_test0.mzml',
                  '../testdata/input/mbr_test1.mzml',
                  '../testdata/input/mbr_test2.mzml')
  
  expected <- arrow::read_parquet('../testdata/unsupervised_recovered_feature_sample_table.parquet')
  
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
  
  result <- unsupervised(test_files, cluster = num_workers)

  report <- dataCompareR::rCompare(result$recovered_feature_sample_table, keys = c("feature", "mz", "rt", "sample_rt", "sample_intensity"), expected, mismatches = 100000)
  dataCompareR::saveReport(report, reportName = "unsupervised_report", showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

  expect_equal(result$recovered_feature_sample_table, expected)
})
