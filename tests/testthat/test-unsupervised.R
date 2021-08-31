test_that("basic unsupervised test", {
  test_files <- c('../testdata/mbr_test0.mzml',
                  '../testdata/mbr_test1.mzml',
                  '../testdata/mbr_test2.mzml')
  
  expected <- arrow::read_parquet('../testdata/unsupervised_recovered_feature_sample_table.parquet')
  
  result <- unsupervised(test_files)

  expect_equal(result$recovered_feature_sample_table, expected)
})
