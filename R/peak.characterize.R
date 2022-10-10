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

characterize <- function(existing.row, n, m, rt.row, ftrs.row){
    existing.row[7] <- sum(existing.row[7], length(ftrs.row) - 4, na.rm = T)
    existing.row[8] <- (n + m) / existing.row[7]
    existing.row[9] <- min(existing.row[6], existing.row[9], ftrs.row[3], na.rm = T)
    existing.row[10] <- max(existing.row[6], existing.row[10], ftrs.row[4], na.rm = T)

    this <- merge.new(existing.row[11], existing.row[12], existing.row[13], existing.row[14], n, rt.row[5:length(rt.row)])
    existing.row[11:14] <- this

    this <- merge.new(existing.row[15], existing.row[16], existing.row[17], existing.row[18], n, ftrs.row[5:length(ftrs.row)])
    existing.row[15:18] <- this

    return(existing.row)
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
#' peak.characterize(existing_row = NA, ftrs_row, rt_row)
peak.characterize <- function(existing_row = NA, ftrs_row, rt_row) {
    ftrs_row[5:length(ftrs_row)] <- log10(ftrs_row[5:length(ftrs_row)] + 1)
    ftrs_row[ftrs_row == 0] <- NA
    if (length(existing_row) == 1) {
        existing_row <- rep(NA, 18)
        existing_row[6] <- ftrs_row[1]
    }

    n <- round(as.numeric(existing_row[7]) * as.numeric(existing_row[8])) # times found in previous experiments
    if (is.na(n)) {
        n <- 0
    }
    m <- sum(!is.na(rt_row[5:length(rt_row)])) # times found in current experiment

    existing_row <- characterize(existing_row, n, m, rt_row, ftrs_row)

    return(tibble::as_tibble(existing_row))
}
