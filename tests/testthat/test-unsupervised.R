test_that("basic unsupervised test", {
  test_files <- c(
    "../testdata/input/mbr_test0.mzml",
    "../testdata/input/mbr_test1.mzml",
    "../testdata/input/mbr_test2.mzml"
  )

  expected <- arrow::read_parquet("../testdata/unsupervised_recovered_feature_sample_table.parquet")

  result <- unsupervised(test_files, cluster = get_num_workers())

  expect_equal(result$recovered_feature_sample_table, expected)
})