patrick::with_parameters_test_that(
  "recover weaker signals test",
  {
    testdata <- file.path("..", "testdata")

    ms_files <- lapply(files, function(x) {
      file.path(testdata, "input", paste0(x, ".mzML"))
    })
    ms_files <- unlist(ms_files)

    filenames <- lapply(files, function(x) {
      file.path(testdata, "extracted", paste0(x, ".parquet"))
    })

    extracted <- lapply(filenames, arrow::read_parquet)
    extracted <- lapply(extracted, as.data.frame)

    filenames <- lapply(files, function(x) {
      file.path(testdata, "adjusted", paste0(x, ".parquet"))
    })

    adjusted <- lapply(filenames, arrow::read_parquet)
    adjusted <- lapply(adjusted, as.data.frame)

    aligned <- load_aligned_features(
      file.path(testdata, "aligned", "rt_cross_table.parquet"),
      file.path(testdata, "aligned", "int_cross_table.parquet"),
      file.path(testdata, "aligned", "tolerances.parquet")
    )

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

    clusterExport(cluster, c(
      "recover.weaker", "load.lcms", "find.turn.point",
      "combine.seq.3", "interpol.area", "duplicate.row.remove", "compute_all_times", "load_file"
    ))
    clusterEvalQ(cluster, library("splines"))

    recovered <- lapply(seq_along(ms_files), function(i) {
      recover.weaker(
        filename = ms_files[[i]],
        sample_name = files[i],
        this.f1 = extracted[[i]],
        this.f2 = adjusted[[i]],
        pk.times = aligned$rt_crosstab,
        aligned.ftrs = aligned$int_crosstab,
        orig.tol = mz_tol,
        align.mz.tol = aligned$mz_tolerance,
        align.chr.tol = aligned$rt_tolerance,
        mz.range = recover_mz_range,
        chr.range = recover_chr_range,
        use.observed.range = use_observed_range,
        bandwidth = 0.5,
        min.bw = min_bandwidth,
        max.bw = max_bandwidth,
        recover.min.count = recover_min_count
      )
    })

    feature_table <- aligned$rt_crosstab[, 1:4]
    rt_crosstab <- cbind(feature_table, sapply(recovered, function(x) x$this.times))
    int_crosstab <- cbind(feature_table, sapply(recovered, function(x) x$this.ftrs))

    feature_names <- rownames(feature_table)
    sample_names <- colnames(aligned$rt_crosstab[, -(1:4)])

    recovered_actual <- list(
      extracted_features = lapply(recovered, function(x) x$this.f1),
      corrected_features = lapply(recovered, function(x) x$this.f2),
      rt_crosstab = as_feature_crosstab(feature_names, sample_names, rt_crosstab),
      int_crosstab = as_feature_crosstab(feature_names, sample_names, int_crosstab)
    )

    aligned_feature_sample_table_actual <- create_feature_sample_table(aligned)
    recovered_feature_sample_table_actual <- create_feature_sample_table(recovered_actual)

    aligned_feature_sample_table_expected <- arrow::read_parquet(file.path(testdata, "recovered", "aligned_feature_sample_table.parquet"))
    recovered_feature_sample_table_expected <- arrow::read_parquet(file.path(testdata, "recovered", "recovered_feature_sample_table.parquet"))

    expect_equal(aligned_feature_sample_table_actual, aligned_feature_sample_table_expected)
    expect_equal(recovered_feature_sample_table_actual, recovered_feature_sample_table_expected)

    # create and load final files

    extracted_recovered_actual <- recovered_actual$extracted_features
    corrected_recovered_actual <- recovered_actual$corrected_features

    filenames <- lapply(files, function(x) {
      file.path(testdata, "recovered", "recovered-extracted", paste0(x, ".parquet"))
    })

    extracted_recovered_expected <- lapply(filenames, arrow::read_parquet)

    filenames <- lapply(files, function(x) {
      file.path(testdata, "recovered", "recovered-corrected", paste0(x, ".parquet"))
    })

    corrected_recovered_expected <- lapply(filenames, arrow::read_parquet)

    # preprocess dataframes
    keys <- c("mz", "pos", "sd1", "sd2")

    extracted_recovered_actual <- lapply(extracted_recovered_actual, function(x) {
      as.data.frame(x) |> dplyr::arrange_at(keys)
    })
    corrected_recovered_actual <- lapply(corrected_recovered_actual, function(x) {
      as.data.frame(x) |> dplyr::arrange_at(keys)
    })

    extracted_recovered_expected <- lapply(extracted_recovered_expected, function(x) {
      as.data.frame(x) |> dplyr::arrange_at(keys)
    })
    corrected_recovered_expected <- lapply(corrected_recovered_expected, function(x) {
      as.data.frame(x) |> dplyr::arrange_at(keys)
    })

    # compare files
    for (i in seq_along(files)) {
      # extracted recovered
      actual_extracted_i <- extracted_recovered_actual[[i]]
      expected_extracted_i <- extracted_recovered_expected[[i]]

      report <- dataCompareR::rCompare(actual_extracted_i, expected_extracted_i, keys = keys, roundDigits = 4, mismatches = 100000)
      dataCompareR::saveReport(report, reportName = files[[i]], showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

      expect_true(nrow(report$rowMatching$inboth) >= 0.9 * nrow(expected_extracted_i))

      incommon <- as.numeric(rownames(report$rowMatching$inboth))

      subset_actual <- actual_extracted_i %>% dplyr::slice(incommon)
      subset_expected <- expected_extracted_i %>% dplyr::slice(incommon)

      expect_equal(subset_actual$area, subset_expected$area, tolerance = 0.01 * max(subset_expected$area))

      # corrected recovered
      actual_corrected_i <- corrected_recovered_actual[[i]]
      expected_corrected_i <- corrected_recovered_expected[[i]]

      report <- dataCompareR::rCompare(actual_corrected_i, expected_corrected_i, keys = keys, roundDigits = 4, mismatches = 100000)
      dataCompareR::saveReport(report, reportName = files[[i]], showInViewer = FALSE, HTMLReport = FALSE, mismatchCount = 10000)

      expect_true(nrow(report$rowMatching$inboth) >= 0.9 * nrow(expected_corrected_i))

      incommon <- as.numeric(rownames(report$rowMatching$inboth))

      subset_actual <- actual_corrected_i %>% dplyr::slice(incommon)
      subset_expected <- expected_corrected_i %>% dplyr::slice(incommon)

      expect_equal(subset_actual$area, subset_expected$area, tolerance = 0.01 * max(subset_expected$area))
    }
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      mz_tol = 1e-05,
      recover_mz_range = NA,
      recover_chr_range = NA,
      use_observed_range = TRUE,
      min_bandwidth = NA,
      max_bandwidth = NA,
      recover_min_count = 3
    )
  )
)

files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened")
mz_tol = 1e-05
recover_mz_range = NA
recover_chr_range = NA
use_observed_range = TRUE
min_bandwidth = NA
max_bandwidth = NA
recover_min_count = 3