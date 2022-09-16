patrick::with_parameters_test_that("basic unsupervised test", { 
  test_files <- sapply(files, function(x) file.path("../testdata/input", x))
  
  expected <- arrow::read_parquet(file.path("../testdata/unsupervised", paste0(.test_name, "_unsupervised.parquet")))
  
  result <- unsupervised(test_files, cluster = get_num_workers())
  keys <- c("mz", "rt", "sample", "sample_rt", "sample_intensity")
  actual <- result$recovered_feature_sample_table

  # # This piece of code serves to re-introduce the actual filenames in the actual outputs.
  # # This was needed for comparison in the previous test cases.
  # levels(actual$sample) <- sapply(files, get_sample_name)
  # actual <- actual |> dplyr::arrange_at(keys)
  # expected <- expected |> dplyr::arrange_at(keys)

  # actual_to_compare <- dplyr::select(actual, keys)
  # expected_to_compare <- dplyr::select(expected, keys)

  report <- dataCompareR::rCompare(actual, expected, keys = keys, roundDigits = 3, mismatches = 100000)
  dataCompareR::saveReport(report, reportName = paste0(.test_name,"_unsupervised_report"), showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

  expect_equal(actual, expected, tolerance = 0.01)
}, patrick::cases(
  mbr_test = list(files = c("mbr_test0.mzml", "mbr_test1.mzml", "mbr_test2.mzml")),
  RCX_shortened = list(files = c("RCX_06_shortened.mzML", "RCX_07_shortened.mzML", "RCX_08_shortened.mzML"))
))
