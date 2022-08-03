find.tol <- function(mz_values,
                     mz_max_diff = 1e-4,
                     aver.bin.size = 4000,
                     min.bins = 50,
                     max.bins = 200,
                     do.plot = TRUE) {
    mz_values <- mz_values[order(mz_values)]
    l <- length(mz_values)
    # pairwise m/z difference divided by their average, filtered outside of tolerance limit
    distances <- (mz_values[2:l] - mz_values[1:(l - 1)]) / 
                 ((mz_values[2:l] + mz_values[1:(l - 1)]) / 2)
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
