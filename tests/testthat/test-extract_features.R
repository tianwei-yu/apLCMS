test_that("extract single feature works", {
  filename <- c("../testdata/mbr_test0.mzml")

  extracted_features <- proc.cdf(filename, cache = FALSE)
  actual <- prof.to.features(extracted_features, do.plot = FALSE)

  expected <- as.data.frame(arrow::read_parquet("../testdata/extracted_features/extracted_0.parquet"))
  actual <- as.data.frame(actual)

  actual <- unique(actual)
  expected <- unique(expected)

  keys <- c("mz", "pos", "sd1", "sd2")

  actual <- dplyr::arrange_at(
    actual, keys
  )
  expected <- dplyr::arrange_at(
    expected, keys
  )

  report <- dataCompareR::rCompare(actual, expected, keys = keys, roundDigits = 4, mismatches = 100000)
  dataCompareR::saveReport(report, reportName = "extract_features_0", showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

  expect_true(nrow(report$rowMatching$inboth) >= 0.9 * nrow(expected))

  incommon <- as.numeric(rownames(report$rowMatching$inboth))

  subset_actual <- actual %>% dplyr::slice(incommon)
  subset_expected <- expected %>% dplyr::slice(incommon)

  expect_equal(subset_actual$area, subset_expected$area, tolerance = 0.01 * max(expected$area))
})
