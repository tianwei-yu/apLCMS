test_that("basic hybrid test", {
  test_files <- c('../testdata/mbr_test0.mzml',
                  '../testdata/mbr_test1.mzml',
                  '../testdata/mbr_test2.mzml')
  
  expected <- arrow::read_parquet('../testdata/hybrid_recovered_feature_sample_table.parquet')
  known_table <- arrow::read_parquet('../testdata/known_table.parquet')

  result <- hybrid(test_files, known_table)

  expect_equal(result$recovered_feature_sample_table, expected)
})
