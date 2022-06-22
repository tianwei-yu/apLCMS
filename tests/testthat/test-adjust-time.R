test_that("adjust time test", {
  extracted_feature_files <- c('../testdata/RCX_01_shortened_v2/RCX_01_shortened_v2_features.Rds',
                               '../testdata/RCX_09_shortened_v2/RCX_09_shortened_v2_features.Rds',
                               '../testdata/RCX_16_shortened_v2/RCX_16_shortened_v2_features.Rds')

  corrected_files <- c('../testdata/adjust-time/RCX_01_shortened_v2_features_corrected.Rds',
                       '../testdata/adjust-time/RCX_09_shortened_v2_features_corrected.Rds',
                       '../testdata/adjust-time/RCX_16_shortened_v2_features_corrected.Rds')

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cluster <- 2L
  } else {
    # use all cores in devtools::test()
    cluster = 4
  }

  if (!is(cluster, 'cluster')) {
    cluster <- parallel::makeCluster(cluster)
    on.exit(parallel::stopCluster(cluster))
  }

  # NOTE: side effect (doParallel has no functionality to clean up)
  doParallel::registerDoParallel(cluster)

  extracted <- lapply(extracted_feature_files, readRDS)

  expected <- lapply(corrected_files, readRDS)
  expected <- lapply(expected, as.data.frame)

  corrected <- adjust.time(
    features = extracted,
    mz.tol = NA,
    chr.tol = NA,
    find.tol.max.d = 10 * 1e-05,
    max.align.mz.diff = 0.01,
    do.plot = FALSE
  )

  corrected <- lapply(corrected, as.data.frame)

  expect_equal(corrected, expected)
})