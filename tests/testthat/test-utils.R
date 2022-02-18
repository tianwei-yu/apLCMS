test_that("splitting files into batches", {
    metadata <- read.table("../testdata/two_step_hybrid_info.csv",
        sep = ",",
        header = TRUE)
    filenames <- c("../testdata/mbr_test0.mzml",
        "../testdata/mbr_test1.mzml",
        "../testdata/mbr_test2.mzml",
        "../testdata/mbr_test0_copy.mzml")


    actual <- split_files_into_batches(
        filenames = filenames,
        metadata = metadata)
    expected <- tibble(filename = filenames, batch = c(1, 1, 2, 2))

    expect_equal(actual, expected)
})