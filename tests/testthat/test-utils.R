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