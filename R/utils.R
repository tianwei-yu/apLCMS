get_feature_values <- function(features, rt_colname) {
    mz <- c()
    chr <- c()
    lab <- c()
    for (i in 1:length(features)) {
        features_batch <- dplyr::as_tibble(features[[i]])
        mz <- c(mz, features_batch$mz)
        chr <- c(chr, features_batch[[rt_colname]])
        lab <- c(lab, rep(i, nrow(features_batch)))
    }
    return(list(mz = mz, chr = chr, lab = lab))
}

extract_pattern_colnames <- function(dataframe, pattern) {
    dataframe <- dplyr::select(dataframe, contains(pattern))
    return(colnames(dataframe))
}