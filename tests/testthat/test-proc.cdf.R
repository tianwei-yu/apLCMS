test_that("test proc.cdf test", {
  filename <- c('../testdata/mbr_test0.mzml')
  actual <- proc.cdf(filename, cache = FALSE)
  expected <- readRDS('../testdata/proc_cdf_expected.Rds')
  expect_equal(actual, expected)
}) 