test_that("basic two-step hybrid test", {
  skip("Disabled")
  skip_on_ci()
  test_names <- c(
    "mbr_test0.mzml",
    "mbr_test1.mzml",
    "mbr_test2.mzml",
    "mbr_test0_copy.mzml"
  )
  test_path <- paste0("../testdata/", test_names)
  metadata <- read.table("../testdata/two_step_hybrid_info.csv", sep = ",", header = TRUE)

  tempdir <- tempdir()
  dir.create(tempdir)
  temp_path <- paste0(tempdir, "/", test_names)
  file.copy(test_path, temp_path)

  expected_final_features <- readRDS("../testdata/final_ftrs.Rda")
  known_table <- arrow::read_parquet("../testdata/known_table.parquet")

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    num_workers <- 2L
  } else {
    num_workers <- parallel::detectCores()
  }

  result <- two.step.hybrid(
    filenames = test_names,
    metadata = metadata,
    work_dir = tempdir,
    known.table = known_table,
    cluster = num_workers
  )
  final_features <- result$final_features

  keys <- c("feature", "mz", "rt", "mz_min", "mz_max", "sample")
  final_features <- as_tibble(arrange_at(final_features, keys))
  expected_final_features <- as_tibble(arrange_at(expected_final_features, keys))
  comparison <- dataCompareR::rCompare(
    final_features,
    expected_final_features,
    keys = keys
  )

  dataCompareR::saveReport(
    comparison,
    reportName = "final_features_comparison",
    reportLocation = ".",
    showInViewer = FALSE,
    missmatchCount = 10000
  )

  unlink(tempdir, recursive = TRUE)

  expect_equal(final_features, expected_final_features, tolerance = 0.001)
  expect_equal(final_features$mz, expected_final_features$mz, tolerance = 0.001)
})
