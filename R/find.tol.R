#' @import MASS
NULL
#> NULL

#' An internal function that is not supposed to be directly accessed by the user. Find m/z tolerance level.
#' 
#' The function finds the tolerance level in m/z from a given vector of observed m/z values.
#' 
#' @param mz_values The vector of observed m/z values.
#' @param mz_max_diff Consider only m/z diffs smaller than this value.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins the maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#' @param do.plot Indicates whether plot should be drawn.
#' @return The tolerance level is returned.
#' @export
#' @examples
#' find.tol(mz_values, mz_max_diff = mz_max_diff, do.plot = FALSE)
find.tol <- function(mz_values,
                     mz_max_diff = 1e-4,
                     aver.bin.size = 4000,
                     min.bins = 50,
                     max.bins = 200,
                     do.plot = TRUE) {
    mz_values <- mz_values[order(mz_values)]
    l <- length(mz_values)
    da <- (mz_values[2:l] - mz_values[1:(l - 1)]) / ((mz_values[2:l] + mz_values[1:(l - 1)]) / 2)
    da <- da[da < mz_max_diff]
    n <- min(max.bins, max(round(length(da) / aver.bin.size), min.bins))
    des <- density(da, kernel = "gaussian", n = n, bw = mz_max_diff / n * 2, from = 0)
    y <- des$y[des$x > 0]
    x <- des$x[des$x > 0]
    
    to.use <- da[da > max(da) / 4] - max(da) / 4
    this.rate <- MASS::fitdistr(to.use, "exponential")$estimate
    exp.y <- dexp(x, rate = this.rate)
    exp.y <- exp.y * sum(y[x > max(da) / 4]) / sum(exp.y[x > max(da) / 4])
    
    yy <- cumsum(y > 1.5 * exp.y)
    yi <- seq_along(yy)
    sel <- min(which(yy < yi)) - 1
    
    if (do.plot) {
        tolerance_plot(x, y, exp.y, sel, main = "find m/z tolerance")
    }
    
    return(x[sel])
}
