patrick::with_parameters_test_that(
  "adjust time test",
  {
    skip_on_ci()
    
    testdata <- file.path("..", "testdata")
    
    filenames <- lapply(files, function(x) {
      file.path(testdata, x, paste0(x, "_features.Rds"))
    })
    
    extracted <- lapply(filenames, readRDS)
    
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
      features = extracted,
      mz_tol_relative = mz_tol,
      rt_tol_relative = chr_tol,
      mz_max_diff = find_tol_max_d,
      mz_tol_absolute = max_align_mz_diff,
      do.plot = FALSE
    )
    
    expected_filenames <- lapply(files, function(x) {
      file.path(testdata, "adjust-time", paste0(x, "_features_corrected.parquet"))
    })
    
    expected <- lapply(expected_filenames, arrow::read_parquet)
    expected <- lapply(expected, as.data.frame)
    
    corrected <- lapply(corrected, as.data.frame)
  
    expect_equal(corrected, expected)
  },
  patrick::cases(
    RCX_shortened = list(
      files = c("RCX_01_shortened_v2", "RCX_09_shortened_v2", "RCX_16_shortened_v2"),
      mz_tol = NA,
      chr_tol = NA,
      find_tol_max_d = 10 * 1e-05,
      max_align_mz_diff = 0.01
    )
  )
)