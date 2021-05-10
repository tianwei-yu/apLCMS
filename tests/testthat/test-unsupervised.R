test_that("basic unsupervised test", {
  test_files <- c('../testdata/mbr_test0.mzml',
                  '../testdata/mbr_test1.mzml',
                  '../testdata/mbr_test2.mzml')

  invisible(capture.output({
    test_result <- unsupervised(test_files)
  }))

  expect_type(test_result, 'list')
})
