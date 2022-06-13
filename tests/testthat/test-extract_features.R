test_that("feature extraction test", {
  ms_files <- c('../testdata/mbr_test0.mzml',
                '../testdata/mbr_test1.mzml',
                '../testdata/mbr_test2.mzml')

  extracted_feature_files <- c('../testdata/extracted_features/extracted_0.parquet',
                               '../testdata/extracted_features/extracted_1.parquet',
                               '../testdata/extracted_features/extracted_2.parquet')

  expected <- lapply(extracted_feature_files, arrow::read_parquet)
  expected <- lapply(expected, as.matrix)

  cluster = 4

  if (!is(cluster, 'cluster')) {
    cluster <- parallel::makeCluster(cluster)
    on.exit(parallel::stopCluster(cluster))
  }

  # NOTE: side effect (doParallel has no functionality to clean up)
  doParallel::registerDoParallel(cluster)

  extracted <- extract_features(
    cluster = cluster,
    filenames = ms_files,
    min_pres = 0.5,
    min_run = 12,
    mz_tol = 1e-05,
    baseline_correct = 0,
    baseline_correct_noise_percentile = 0.05,
    intensity_weighted = FALSE,
    min_bandwidth = NA,
    max_bandwidth = NA,
    sd_cut = c(0.01, 500),
    sigma_ratio_lim = c(0.01, 100),
    shape_model = "bi-Gaussian",
    peak_estim_method = "moment",
    component_eliminate = 0.01,
    moment_power = 1,
    BIC_factor = 2
  )

  expect_equal(extracted, expected)
})
