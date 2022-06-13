test_that("prof.to.features works", {
    extracted_features <- readRDS('../testdata/proc_cdf_expected.Rds')
    actual <- prof.to.features(extracted_features, do.plot = FALSE)

    expected <- readRDS("../testdata/prof_to_features_expected.Rds")
    expect_equal(actual, expected)
})