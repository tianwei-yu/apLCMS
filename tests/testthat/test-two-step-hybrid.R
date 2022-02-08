test_that("basic two-step hybrid test", {
    test_folder <- "../testdata"

    expected <- "..."
    known_table <- arrow::read_parquet('../testdata/known_table.parquet')
    info <- read.table("../testdata/two_step_hybrid_info.csv", sep=",")

    # CRAN limits the number of cores available to packages to 2
    # source https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions#50571533
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
    if (nzchar(chk) && chk == "TRUE") {
        # use 2 cores in CRAN/Travis/AppVeyor
        num_workers <- 2L
    } else {
        # use all cores in devtools::test()
        num_workers <- parallel::detectCores()
    }

    result <- two.step.hybrid(folder=test_folder, info=info, file.pattern="mzml", known.table=known_table)
    saveRDS(result, "result.Rda")

    expect_equal(result$recovered_feature_sample_table, expected)

})