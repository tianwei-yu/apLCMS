test_that("basic hybrid test", {
  test_files <- c('../testdata/mbr_test0.mzml',
                  '../testdata/mbr_test1.mzml',
                  '../testdata/mbr_test2.mzml')

  capture.output({
    test_result <- hybrid(test_files, known.table.hplus)
  }, type = 'message')

  expect_type(test_result, 'list')
})
