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
