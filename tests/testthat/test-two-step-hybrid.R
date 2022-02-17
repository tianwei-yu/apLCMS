test_that("basic two-step hybrid test", {
    test_files <- c('../testdata/mbr_test0.mzml',
                    '../testdata/mbr_test1.mzml',
                    '../testdata/mbr_test2.mzml',
                    '../testdata/mbr_test0_copy.mzml')
    metadata <- "../testdata/two_step_hybrid_info.csv"

    expected_final_ftrs <- readRDS("../testdata/final_ftrs.Rda")
    expected_all_ftrs <- readRDS("../testdata/all.detected.ftrs.Rda")
    expected_batchwise <- readRDS("../testdata/batchwise.results.Rda")
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

    result <- two.step.hybrid(filenames=test_files, metadata=metadata, known.table=known_table)

    expect_equal(result$final.ftrs, expected_final_ftrs)
    expect_equal(result$all.detected.ftrs, expected_all_ftrs)
    expect_equal(result$batchwise.results, expected_batchwise)

})