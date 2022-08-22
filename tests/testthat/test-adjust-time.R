patrick::with_parameters_test_that(
  "adjust time test",
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
      mz_tol_relative = mz_tol_relative,
      rt_tol_relative = rt_tol_relative,
      mz_max_diff = 10 * mz_tol,
      mz_tol_absolute = mz_tol_absolute,
      do.plot = do.plot 
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
      mz_tol_relative = NA,
      rt_tol_relative = NA,
      mz_tol = 1e-05,
      mz_tol_absolute = 0.01,
      do.plot = FALSE
    )
  )
)
