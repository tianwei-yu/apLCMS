test_that("adjust time test", {
  extracted_feature_files <- c('../testdata/extracted-feature/extracted_0.parquet',
                               '../testdata/extracted-feature/extracted_1.parquet',
                               '../testdata/extracted-feature/extracted_2.parquet')
  
  corrected_files <- c('../testdata/adjust-time/corrected_0.parquet',
                       '../testdata/adjust-time/corrected_1.parquet',
                       '../testdata/adjust-time/corrected_2.parquet')

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cluster <- 2L
  } else {
    # use all cores in devtools::test()
    cluster <- parallel::detectCores()
  }

  if (!is(cluster, 'cluster')) {
    cluster <- parallel::makeCluster(cluster)
    on.exit(parallel::stopCluster(cluster))
  }

  # NOTE: side effect (doParallel has no functionality to clean up)
  doParallel::registerDoParallel(cluster)

  extracted <- lapply(extracted_feature_files, arrow::read_parquet)
  extracted <- lapply(extracted, as.matrix)

  expected <- lapply(corrected_files, arrow::read_parquet)
  expected <- lapply(expected, as.matrix)

  corrected <- adjust.time(
    features = extracted,
    mz.tol = NA,
    chr.tol = NA,
    find.tol.max.d = 10 * 1e-05,
    max.align.mz.diff = 0.01,
    do.plot = FALSE
  )

  expect_equal(corrected, expected)
})