#' An internal function.
#'
#' This is a internal function.
#'
#' @param features dataframe of retention time, m/z ratio, signal strength.
#' @return returns
#' \itemize{
#'   \item masses - m/z ratio
#'   \item labels - retention time
#'   \item intensi - signal strength
#' }
#' @export
#' @examples
#' combine.seq.3(table)
combine.seq.3 <- function(features) {
    l <- nrow(features)
    breaks <- compute_breaks_3(features$labels)
    new_table <- tibble::tibble(
        mz = rep(0, length(breaks) - 1),
        labels = unique(features$labels),
        intensities = rep(0, length(breaks) - 1)
    )

    for (i in 1:(length(breaks) - 1)) {
        start <- breaks[i] + 1
        end <- breaks[i + 1]

        mz <- features$mz[start:end]
        ints <- features$intensities[start:end]

        new_table$intensities[i] <- sum(ints)
        new_table$mz[i] <- median(mz[which.max(ints)])
    }

    return(new_table)
}
