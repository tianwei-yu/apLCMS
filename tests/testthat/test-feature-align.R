test_that("feature align test", {
  ms_files <- c('../testdata/mbr_test0.mzml',
                '../testdata/mbr_test1.mzml',
                '../testdata/mbr_test2.mzml')
  files_list <- sort_samples_by_acquisition_number(ms_files)
  sample_names <- get_sample_name(files_list)
  
  corrected_files <- c('../testdata/feature-align/corrected_0.parquet',
                       '../testdata/feature-align/corrected_1.parquet',
                       '../testdata/feature-align/corrected_2.parquet')
  
  files_list <- sort_samples_by_acquisition_number(corrected_files)
  corrected_features <- lapply(files_list, arrow::read_parquet)
  corrected_features <- lapply(corrected_features, as.matrix)
  
  aligned_actual <- align_features(
      sample_names = sample_names,
      features = corrected_features,
      min.exp = 2,
      mz.tol = NA,
      chr.tol = NA,
      find.tol.max.d = 10 * 1e-05,
      max.align.mz.diff = 0.01,
      do.plot = FALSE
  )
  
  rt_cross_table <- arrow::read_parquet('../testdata/feature-align/rt_cross_table.parquet')
  int_cross_table <- arrow::read_parquet('../testdata/feature-align/int_cross_table.parquet')
  tolerances_table <- arrow::read_parquet('../testdata/feature-align/tolerances.parquet')
  
  row.names(rt_cross_table) <- as.character(row.names(rt_cross_table))
  row.names(int_cross_table) <- as.character(row.names(int_cross_table))
  
  aligned_expected <- list()
  aligned_expected$mz_tolerance <- tolerances_table$mz_tolerance
  aligned_expected$rt_tolerance <- tolerances_table$rt_tolerance
  aligned_expected$rt_crosstab <- rt_cross_table
  aligned_expected$int_crosstab <- int_cross_table

  expect_equal(aligned_actual,aligned_expected)
})