#' @import MASS
NULL
#> NULL

#' An internal function that is not supposed to be directly accessed by the user. Find m/z tolerance level.
#' 
#' The function finds the tolerance level in m/z from a given vector of observed m/z values.
#' 
#' @param mz The vector of observed m/z values.
#' @param mz_max_diff Consider only m/z diffs smaller than this value.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins the maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#' @param do.plot Indicates whether plot should be drawn.
#' @return The tolerance level is returned.
#' @export
#' @examples
#' find.tol(mz, mz_max_diff = mz_max_diff, do.plot = FALSE)
find.tol <- function(mz,
                     mz_max_diff = 1e-4,
                     aver.bin.size = 4000,
                     min.bins = 50,
                     max.bins = 200,
                     do.plot = TRUE) {
    mz <- mz[order(mz)]
    l <- length(mz)
    # pairwise m/z difference divided by their average, filtered outside of tolerance limit
    distances <- (mz[2:l] - mz[1:(l - 1)]) / 
                 ((mz[2:l] + mz[1:(l - 1)]) / 2)
    distances <- distances[distances < mz_max_diff]
    
    # number of equally spaced points at which the density is to be estimated
    n <- min(max.bins, max(round(length(distances) / aver.bin.size), min.bins))
    # estimate probability density function of distances
    des <- density(distances, kernel = "gaussian", n = n,
                   bw = mz_max_diff / n * 2, from = 0)
    # the n (-1?) coordinates of the points where the density is estimated
    points <- des$y[des$x > 0]
    # the estimated density values
    density_values <- des$x[des$x > 0]
    
    # select the upper 75% of the sorted data
    top_data <- distances[distances > max(distances) / 4] - max(distances) / 4
    # parameter of the exponential distribution is estimated
    lambda <- MASS::fitdistr(top_data, "exponential")$estimate
    # values of the exponential distribution are calculated at equally spaced points
    estimated_density_values <- dexp(density_values, rate = lambda)
    # add the rest of the data
    estimated_density_values <- estimated_density_values * sum(points[density_values > max(distances) / 4]) / 
                                                           sum(estimated_density_values[density_values > max(distances) / 4])
    
    # cutoff is selected where the density of the empirical distribution is >1.5 times the density of the exponential distribution
    cumulative <- cumsum(points > 1.5 * estimated_density_values)
    cumulative_indices <- seq_along(cumulative)
    selected <- min(which(cumulative < cumulative_indices)) - 1
    
    if (do.plot) {
        tolerance_plot(density_values, points, estimated_density_values, selected, main = "find m/z tolerance")
    }
    
    return(density_values[selected])
}
