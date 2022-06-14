test_that("prof.to.features works", {
    extracted_features <- readRDS('../testdata/proc_cdf_expected.Rds')
    actual <- prof.to.features(
        extracted_features,
        sd.cut = c(0.1, 100),
        sigma.ratio.lim = c(0.1, 10),
        do.plot = FALSE)

    expected <- readRDS("../testdata/prof_to_features_expected.Rds")
    expect_equal(actual, expected)
})