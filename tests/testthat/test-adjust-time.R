patrick::with_parameters_test_that(
  "get_template",
  {
    testdata <- file.path("..", "testdata")
    
    filenames <- lapply(files, function(x) {
      file.path(testdata, "clusters", paste0(x, "_extracted_clusters.parquet"))
    })
    
    extracted <- lapply(filenames, arrow::read_parquet)

    template_features <- compute_template(extracted)

    expected <- file.path(testdata, "template", "RCX_shortened.parquet")
    expected <- arrow::read_parquet(expected)

    expect_equal(template_features, expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened")
    )
  )
)

patrick::with_parameters_test_that(
  "adjust time test",
  {
    testdata <- file.path("..", "testdata")

    template_features <- file.path(testdata, "template", "RCX_shortened.parquet")
    template_features <- arrow::read_parquet(template_features)

    extracted <- file.path(testdata, "clusters", paste0(.test_name, "_extracted_clusters.parquet"))
    extracted <- arrow::read_parquet(extracted)

    corrected <- correct_time(
      this.feature = extracted,
      template_features = template_features,
      mz_tol_relative = mz_tol_relative,
      rt_tol_relative = rt_tol_relative
    )

    expected <- tibble::as_tibble(arrow::read_parquet(file.path(testdata, "adjusted", paste0(.test_name, ".parquet"))))
    expected <- expected |> dplyr::rename( rt = pos, sample_id = V6, cluster = V7)

    expect_equal(corrected, expected)
  },
  patrick::cases(
    RCX_06_shortened = list(
      mz_tol_relative = 6.856763e-06,
      rt_tol_relative = 3.618581
    ),
    RCX_07_shortened = list(
      mz_tol_relative = 6.856763e-06,
      rt_tol_relative = 3.618581
    ),
    RCX_08_shortened = list(
      mz_tol_relative = 6.856763e-06,
      rt_tol_relative = 3.618581
    )
  )
)


patrick::with_parameters_test_that(
  "adjust time complete",
  {
    testdata <- file.path("..", "testdata")
    
    filenames <- lapply(files, function(x) {
      file.path(testdata, "extracted", paste0(x, ".parquet"))
    })
    
    extracted <- lapply(filenames, arrow::read_parquet)
    extracted <- lapply(extracted, as.data.frame)
    
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
    
    # NOTE: side effect (doParallel has no functionality to clean up)
    doParallel::registerDoParallel(cluster)
    
    corrected <- adjust.time(
      extracted_features = extracted,
      mz_tol_relative = mz_tol,
      rt_tol_relative = chr_tol,
      mz_max_diff = find_tol_max_d,
      mz_tol_absolute = max_align_mz_diff,
      do.plot = FALSE
    )
    
    expected_filenames <- lapply(files, function(x) {
      file.path(testdata, "adjusted", paste0(x, ".parquet"))
    })
    
    expected <- lapply(expected_filenames, function(x) {
      tibble::as_tibble(arrow::read_parquet(x)) |> dplyr::rename( rt = pos, sample_id = V6, cluster = V7)
    })
    
    corrected <- lapply(corrected, tibble::as_tibble)
  
    expect_equal(corrected, expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_06_shortened", "RCX_07_shortened", "RCX_08_shortened"),
      mz_tol = NA,
      chr_tol = NA,
      find_tol_max_d = 10 * 1e-05,
      max_align_mz_diff = 0.01
    )
  )
)