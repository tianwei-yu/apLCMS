patrick::with_parameters_test_that(
  "basic unsupervised test",
  {
    store_reports <- FALSE

    if (skip) {
      skip("skipping whole data test case")
    }

    test_files <- sapply(files, function(x) file.path("../testdata/input", x))

    # expected <- arrow::read_parquet(file.path("../testdata/unsupervised", paste0(.test_name, "_unsupervised.parquet")))

    result <- unsupervised(test_files, cluster = get_num_workers())
    keys <- c("mz", "rt", "sample", "sample_rt", "sample_intensity")
    actual <- result$recovered_feature_sample_table

    if (store_reports) {
      report <- dataCompareR::rCompare(
        actual,
        expected,
        keys = keys,
        roundDigits = 3,
        mismatches = 100000
      )
      dataCompareR::saveReport(
        report,
        reportName = paste0(.test_name, "_unsupervised_report"),
        showInViewer = FALSE,
        HTMLReport = FALSE,
        mismatchCount = 10000
      )
    }

    # write_parquet(actual, "qc_no_dil_milliq_unsupervised.parquet")

    # expect_equal(actual, expected, tolerance = 0.01)
  },
  patrick::cases(
    mbr_test = list(
      files = c("mbr_test0.mzml", "mbr_test1.mzml", "mbr_test2.mzml"),
      skip = TRUE
    ),
    RCX_shortened = list(
      files = c("RCX_06_shortened.mzML", "RCX_07_shortened.mzML", "RCX_08_shortened.mzML"),
      skip = TRUE
    ),
    qc_no_dil_milliq = list(
      files = c("8_qc_no_dil_milliq.mzml", "21_qc_no_dil_milliq.mzml", "29_qc_no_dil_milliq.mzml"),
      skip = FALSE
    )
  )
)
