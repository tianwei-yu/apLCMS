patrick::with_parameters_test_that(
  "test benchmark",
  {
    if (skip) {
      skip("Disabled")
    }

    skip_on_ci()

    testdata <- file.path("..", "testdata")

    filenames <- lapply(filename, function(x) {
      file.path(testdata, "input", paste0(x, ".mzml"))
    })

    cluster <- get_num_workers()

    if (!is(cluster, "cluster")) {
      cluster <- parallel::makeCluster(cluster)
      on.exit(parallel::stopCluster(cluster))
    }

    res <- microbenchmark::microbenchmark(
      unsupervised = {
        result <- unsupervised(unlist(filenames), cluster = cluster)
      },
      times = 10L
    )

    expected_path <- file.path(testdata, expected_filename)
    expected <- arrow::read_parquet(expected_path)
    expect_equal(result$recovered_feature_sample_table, expected)

    cat("\n\n")
    print(res)
    cat("\n")
  },
  patrick::cases(
    RCX_shortened = list(
      filename = c("mbr_test0", "mbr_test1", "mbr_test2"),
      expected_filename = "unsupervised_recovered_feature_sample_table.parquet",
      skip = TRUE
    )
  )
)
