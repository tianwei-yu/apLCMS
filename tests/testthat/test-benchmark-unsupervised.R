patrick::with_parameters_test_that(
  "test benchmark",
  {
    if(skip){
      skip("Disabled")
    }

    skip_on_ci()

    testdata <- file.path("..", "testdata")

    filenames <- lapply(filename, function(x) {
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
