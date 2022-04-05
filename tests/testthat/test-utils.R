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

test_that("Utils: colnames pattern extraction", {
    dataframe <- dplyr::tibble(
        mz = numeric(),
        rt = numeric(),
        sample_1_intensity = numeric(),
        sample_2_intensity = numeric(),
        sample_1_rt = numeric(),
        sample_2_rt = numeric()
    )

    intensity_labels <- extract_pattern_colnames(dataframe, "_intensity")
    rt_labels <- extract_pattern_colnames(dataframe, "_rt")

    expect_equal(intensity_labels, c("sample_1_intensity", "sample_2_intensity"))
    expect_equal(rt_labels, c("sample_1_rt", "sample_2_rt"))
})