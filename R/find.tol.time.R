#' An internal function that find elution time tolerance level. 
#' 
#' This function finds the time tolerance level. Also, it returns the grouping information given the time tolerance.
#' 
#' @param mz mz values of all peaks in all profiles in the study.
#' @param chr retention time of all peaks in all profiles in the study.
#' @param lab label of all peaks in all profiles in the study.
#' @param number_of_samples The number of spectra in this analysis.
#' @param mz_tol_relative m/z tolerance level for the grouping of signals into peaks. This value is expressed as the percentage of the m/z value. 
#'  This value, multiplied by the m/z value, becomes the cutoff level.
#' @param rt_tol_relative the elution time tolerance. If NA, the function finds the tolerance level first. If a numerical value is given, 
#'  the function directly goes to the second step - grouping peaks based on the tolerance.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins the maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#' @param mz_tol_absolute As the m/z tolerance in alignment is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. 
#'  This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.
#' @param max.num.segments the maximum number of segments.
#' @param do.plot Indicates whether plot should be drawn.
#' @return A matrix with six columns. Every row corrsponds to a peak in one of the spectrum. The columns are: m/z, elution time, spread, signal strength, 
#'  spectrum label, and peak group label. The rows are ordered by the median m/z of each peak group, and with each peak group the rows are ordered 
#'  by the elution time.
#' @examples
#' find.tol.time(mz, chr, lab, number_of_samples = number_of_samples, mz_tol_relative = mz_tol_relative, mz_tol_absolute = mz_tol_absolute, do.plot = FALSE)
find.tol.time <- function(mz,
                          rt,
                          sample_id,
                          number_of_samples,
                          mz_tol_relative = 2e-5,
                          rt_tol_relative = NA,
                          aver.bin.size = 200,
                          min.bins = 50,
                          max.bins = 100,
                          mz_tol_absolute = 0.01,
                          max.num.segments = 10000,
                          do.plot = TRUE) {
    # sort m/z, rt, and sampel IDs based on m/z
    mz_ordering <- order(mz)
    mz <- mz[mz_ordering]
    rt <- rt[mz_ordering]
    sample_id <- sample_id[mz_ordering]
    rm(mz_ordering)
    
    l <- length(mz)
    
    # indices of m/z where difference from its neighbor is greater than allowed
    # tolerance; as result, it separates m/z to groups with similar values
    breaks <-
      c(0, which((mz[2:l] - mz[1:(l - 1)]) 
                  > min(mz_tol_absolute, 
                       mz_tol_relative * ((mz[2:l] + mz[1:(l - 1)]) / 2)
                       )
                ), 
        l)
    
    for (i in 2:length(breaks)) {
        # sort rt inside each m/z group
        rt_ordering <- order(rt[(breaks[i - 1] + 1):breaks[i]])
        rt_ordering <- rt_ordering + breaks[i - 1]
        # reorder other m/z, rt and sample ID within group based on rt order
        mz[(breaks[i - 1] + 1):breaks[i]] <- mz[rt_ordering]
        rt[(breaks[i - 1] + 1):breaks[i]] <- rt[rt_ordering]
        sample_id[(breaks[i - 1] + 1):breaks[i]] <- sample_id[rt_ordering]
    }
    
    # remove fist and last index
    breaks <- breaks[c(-1, -length(breaks))]
    
    if (is.na(rt_tol_relative)) {
        distances <- 0
        # create indices for each groups
        if (length(breaks) > max.num.segments) {
            segments <- floor(seq(2, length(breaks), length.out = max.num.segments))
        } else {
            segments <- 2:length(breaks)
        }
        
        # compute distance matrix of each group using stats::dist
        for (i in segments) {
            segment <- (breaks[i - 1] + 1):breaks[i]
            
            if (length(segment) <= 3 * number_of_samples) {
                distance_matrix <- as.vector(dist(rt[segment]))
                if (length(distance_matrix) > 100)
                    distance_matrix <- sample(distance_matrix, 100)
                distances <- c(distances, distance_matrix)
            }
        }
        
        # goal: statistically find the smallest rt which is still significant
        
        # a long vector of distances between rt values (with no particular order)
        distances <- distances[!is.na(distances)]
        max_distance <- max(distances) # maximal distance
        # number of equally spaced points at which the density is to be estimated
        n = min(max.bins, max(min.bins, round(length(distances) / aver.bin.size)))
        # estimate probability density function of distances
        des <- density(distances, kernel = "gaussian", n = n,
                       bw = max_distance / n * 2, from = 0)
        # the n (-1?) coordinates of the points where the density is estimated
        points <- des$x[des$x > 0]
        # the estimated density values
        density_values <- des$y[des$x > 0]
        
        linear_regression_model <- lm(density_values[points > max_distance / 4] ~ points[points > max_distance / 4])
        
        # compute probability density values (y) using the linear function
        estimated_density_values <- linear_regression_model$coef[1] + linear_regression_model$coef[2] * points
        
        values_not_last <- density_values[1:(length(density_values) - 1)] # density values without the last one
        values_not_first <- density_values[2:(length(density_values))]    # density values without the first one
        # pair-wise copy greater value to the left
        values_not_last[which(values_not_last < values_not_first)] <- values_not_first[which(values_not_last < values_not_first)]
        density_values[1:(length(density_values) - 1)] <- values_not_last
        
        # cumulative sum - where density value is greater than estimated density value
        # cutoff is selected where the density of the empirical distribution is >1.5 times the density of the distribution
        cumulative <- cumsum(density_values > 1.5 * estimated_density_values)
        cumulative_indices <- seq_along(cumulative)
        # find last index where density value is greater than estimated density value
        selected <- min(which(cumulative < cumulative_indices)) - 1
        # corresponding coordinate is used as rt tolerance
        rt_tol_relative <- points[selected]
        
        if (do.plot) {
          tolerance_plot(points, density_values, estimated_density_values, selected, main = "find retention time tolerance")
        }
    }
    
    # pair-wise rt differences
    distances <- rt[2:l] - rt[1:(l - 1)]
    # find those respecting the computed tolerance
    breaks_rt <- which(distances > rt_tol_relative)
    # merge and sort both group delimiters 
    breaks_final <- c(0, unique(c(breaks, breaks_rt)), l)
    breaks_final <- breaks_final[order(breaks_final)]
    
    # create list of indices corresponding to a representative from each group 
    # (always the first element)
    groups <- seq_along(mz)
    for (i in 2:length(breaks_final)) {
        groups[(breaks_final[i - 1] + 1):breaks_final[i]] <- i
    }
    
    list(
        mz = mz,
        rt = rt,
        lab = sample_id,
        grps = groups,
        rt.tol = rt_tol_relative,
        mz.tol = mz_tol_relative
    )
}
