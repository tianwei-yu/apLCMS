patrick::with_parameters_test_that("basic unsupervised test", { 
  test_files <- sapply(files, function(x) file.path("../testdata/input", x))
  
  expected <- arrow::read_parquet(file.path("../testdata/recovered", paste0(.test_name, ".parquet")))
  
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
  keys <- c("mz", "rt", "sample", "sample_rt", "sample_intensity")
  actual <- result$recovered_feature_sample_table

  levels(actual$sample) <- sapply(files, get_sample_name)
  actual <- actual |> dplyr::arrange_at(keys)
  expected <- expected |> dplyr::arrange_at(keys)

  actual_to_compare <- dplyr::select(actual, keys)
  expected_to_compare <- dplyr::select(expected, keys)

  report <- dataCompareR::rCompare(actual, expected, keys = keys, roundDigits = 3, mismatches = 100000)
  dataCompareR::saveReport(report, reportName = paste0(.test_name,"_unsupervised_report"), showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

  expect_equal(actual_to_compare, expected_to_compare, tolerance = 0.01)
}, patrick::cases(
  mbr_test = list(files = c("mbr_test0.mzml", "mbr_test1.mzml", "mbr_test2.mzml")),
  RCX_shortened = list(files = c("RCX_06_shortened.mzML", "RCX_07_shortened.mzML", "RCX_08_shortened.mzML"))
))
