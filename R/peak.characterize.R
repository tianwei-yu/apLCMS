merge.new <- function(mean0, sd0, min0, max0, n, x) {
    x <- x[!is.na(x)]
    if (n <= 1) {
        if (n == 1) {
            x <- c(x, mean0)
        }
        mean1 <- mean(x)
        sd1 <- sd(x)
        min1 <- min(x, Inf)
        max1 <- max(x, -Inf)
    } else {
        m <- length(x)
        mean1 <- sum(mean0 * n, x) / (n + m)
        sd1 <- sqrt((n * (mean0 - mean1)^2 + sum((x - mean1)^2) + (n - 1) * sd0^2) / (m + n - 1))
        min1 <- min(min0, x)
        max1 <- max(max0, x)
    }
    return(c(mean1, sd1, min1, max1))
}

characterize <- function(existing_row, n, m, metadata_row, rt_row, ftrs_row) {
    existing_row[7] <- sum(existing_row[7], length(ftrs_row) - 1, na.rm = T)
    existing_row[8] <- (n + m) / existing_row[7]
    existing_row[9] <- min(existing_row[6], existing_row[9], metadata_row$mzmin, na.rm = T)
    existing_row[10] <- max(existing_row[6], existing_row[10], metadata_row$mzmax, na.rm = T)

    this <- merge.new(existing_row[11], existing_row[12], existing_row[13], existing_row[14], n, rt_row[2:length(rt_row)])
    existing_row[11:14] <- this

    this <- merge.new(existing_row[15], existing_row[16], existing_row[17], existing_row[18], n, ftrs_row[2:length(ftrs_row)])
    existing_row[15:18] <- this

    return(existing_row)
}

#' Internal function: Updates the information of a feature for the known feature table.
#'
#' @description
#' The function takes the information about the feature in the known feature table (if available), and updates it using the
#' information found in the current dataset.
#' @param existing_row The existing row in the known feature table.
#' @param ftrs_row The row of the matched feature in the new aligned feature table.
#' @param rt_row The row of the matched feature in the new retention time table of aligned features.
#' @return A vector, the updated row for the known feature table.
#' @examples
#' peak_characterize(existing_row = NA, ftrs_row, rt_row)
peak_characterize <- function(existing_row = NA, metadata_row, ftrs_row, rt_row) {
    ftrs_row[2:length(ftrs_row)] <- log10(ftrs_row[2:length(ftrs_row)] + 1)
    ftrs_row[ftrs_row == 0] <- NA
    if (length(existing_row) == 1) {
        existing_row <- rep(NA, 18)
        existing_row[6] <- metadata_row$mz
    }

    n <- round(as.numeric(existing_row[7]) * as.numeric(existing_row[8])) # times found in previous experiments
    if (is.na(n)) {
        n <- 0
    }
    m <- sum(!is.na(rt_row[2:length(rt_row)])) # times found in current experiment

    existing_row <- characterize(existing_row, n, m, metadata_row, rt_row, ftrs_row)

    return(tibble::as_tibble(existing_row))
}
