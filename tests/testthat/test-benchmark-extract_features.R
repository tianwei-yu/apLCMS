patrick::with_parameters_test_that(
  "test benchmark",
  {
    if (skip) {
      skip("Disabled")
    }

    skip_on_ci()

    testdata <- file.path("..", "testdata")

    filenames <- lapply(filename, function(x) {
      file.path(testdata, "input", paste0(x, ".mzML"))
    })

    cluster <- get_num_workers()

    if (!is(cluster, "cluster")) {
      cluster <- parallel::makeCluster(cluster)
      on.exit(parallel::stopCluster(cluster))
    }

    res <- microbenchmark::microbenchmark(
      extract_feature = {
        actual <- extract_features(
          cluster = cluster,
          filenames,
          min_presence = min_presence,
          min_elution_length = min_elution_length,
          mz_tol = mz_tol,
          baseline_correct = 0,
          baseline_correct_noise_percentile = 0.05,
          intensity_weighted = intensity_weighted,
          min_bandwidth = NA,
          max_bandwidth = NA,
          sd_cut = sd_cut,
          sigma_ratio_lim = sigma_ratio_lim,
          shape_model = "bi-Gaussian",
          peak_estim_method = "moment",
          component_eliminate = 0.01,
          moment_power = 1,
          BIC_factor = 2.0
        )
      },
      times = 10L
    )

    expected_filenames <- lapply(filename, function(x) {
      file.path(testdata, "extracted", paste0(x, ".parquet"))
    })
    expected <- lapply(expected_filenames, arrow::read_parquet)
    expected <- lapply(expected, as.data.frame)
    actual <- unique(actual)
    expected <- unique(expected)
    keys <- c("mz", "pos", "sd1", "sd2")
    actual <- lapply(actual, function(x) {
      as.data.frame(x) |> dplyr::arrange_at(keys)
    })
    expected <- lapply(expected, function(x) {
      as.data.frame(x) |> dplyr::arrange_at(keys)
    })
    for (i in seq_along(filename)) {
      actual_i <- actual[[i]]
      expected_i <- expected[[i]]
      report <- dataCompareR::rCompare(actual_i, expected_i, keys = keys, roundDigits = 4, mismatches = 100000)
      dataCompareR::saveReport(report, reportName = filename[[i]], showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)
      expect_true(nrow(report$rowMatching$inboth) >= 0.9 * nrow(expected_i))
      incommon <- as.numeric(rownames(report$rowMatching$inboth))
      subset_actual <- actual_i %>% dplyr::slice(incommon)
      subset_expected <- expected_i %>% dplyr::slice(incommon)
      expect_equal(subset_actual$area, subset_expected$area, tolerance = 0.01 * max(subset_expected$area))
    }

    cat("\n\n")
    print(res)
    cat("\n")
  },
  patrick::cases(
    RCX_shortened = list(
      filename = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      mz_tol = 1e-05,
      min_presence = 0.5,
      min_elution_length = 12,
      intensity_weighted = FALSE,
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      skip = TRUE
    )
  )
)
