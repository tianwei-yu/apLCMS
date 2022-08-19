patrick::with_parameters_test_that(
  "test benchmark",
  {
    skip_on_ci()

    testdata <- file.path("..", "testdata")

    filenames_extract <- lapply(files_extract, function(x) {
      file.path(testdata, "input", paste0(x, ".mzML"))
    })

    filenames_unsupervised <- lapply(files_unsupervised, function(x) {
      file.path(testdata, paste0(x, ".mzml"))
    })

    # CRAN limits the number of cores available to packages to 2
    # source https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions#50571533
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      cluster <- 2L
    } else {
      # use all cores in devtools::test()
      cluster <- parallel::detectCores()
    }

    if (!is(cluster, "cluster")) {
      cluster <- parallel::makeCluster(cluster)
      on.exit(parallel::stopCluster(cluster))
    }

    if(skip){
      skip("Disabled")
    }

    res <- microbenchmark::microbenchmark(
      unsupervised = {
          result <- unsupervised(unlist(filenames_unsupervised), cluster = cluster)
      },
      extract_feature = {
          actual <- extract_features(
            filenames_extract,
            cluster = cluster,
            min_pres = min_pres,
            min_run = min_run,
            mz_tol = tol,
            intensity_weighted = intensity_weighted,
            sd_cut = sd_cut,
            sigma_ratio_lim = sigma_ratio_lim,
            baseline_correct = 0,
            baseline_correct_noise_percentile = 0.05,
            min_bandwidth = NA,
            max_bandwidth = NA,
            moment_power = 1,
            BIC_factor = 2.0,
            component_eliminate = 0.01,
            peak_estim_method = "moment",
            shape_model = "bi-Gaussian"
          )
      },
      times = 10L
    )

    expected <- arrow::read_parquet('../testdata/unsupervised_recovered_feature_sample_table.parquet')
    expect_equal(result$recovered_feature_sample_table, expected)

    expected_extract_filenames <- lapply(files_extract, function(x) {
      file.path(testdata, "extracted", paste0(x, ".parquet"))
    })
    expected <- lapply(expected_extract_filenames, arrow::read_parquet)
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
    for (i in seq_along(files_extract)) {
      actual_i <- actual[[i]]
      expected_i <- expected[[i]]
      report <- dataCompareR::rCompare(actual_i, expected_i, keys = keys, roundDigits = 4, mismatches = 100000)
      dataCompareR::saveReport(report, reportName = files_extract[[i]], showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)
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
      files_extract = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      files_unsupervised = c("mbr_test0", "mbr_test1", "mbr_test2"),
      tol = 1e-05,
      min_pres = 0.5,
      min_run = 12,
      intensity_weighted = FALSE,
      sd_cut = c(0.01, 500),
      sigma_ratio_lim = c(0.01, 100),
      skip = TRUE
    )
  )
)
