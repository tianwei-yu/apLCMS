#' Compute minimum mz tolerance to use.
#' @description
#' Compute the minimum mz tolerance based on the relative
#' tolerance and the mz values and the absolute tolerance.
#' Uses midpoints between mz values for the weighting.
#' @param mz vector Mz values to use.
#' @param mz_tol_relative float Relative mz tolerance to use with the mz values.
#' This forms a sort of weighted tolerance.
#' @param mz_tol_absolute float Absolute tolerance to use independent from the mz values.
#' @return float Minimum tolerance values to use.

compute_min_mz_tolerance <- function(mz, mz_tol_relative, mz_tol_absolute) {
    l <- length(mz)
    mz_midpoints <- ((mz[2:l] + mz[1:(l - 1)]) / 2)
    mz_ftr_relative_tolerances <- mz_tol_relative * mz_midpoints
    min_mz_tol <- min(mz_tol_absolute, mz_ftr_relative_tolerances)
    return(min_mz_tol)
}

#' @description
#' Compute indices of mass differences greater than min_mz_tol.
#' @param mz mz values of all peaks in all profiles in the study.
#' @param min_mz_tol float Minimum tolerance value.
#' @return breaks Integer indices of mass differences to use.
compute_breaks_3 <- function(mz, min_mz_tol) {
    l <- length(mz)
    mass_differences <- diff(mz)
    indices <- which(mass_differences > min_mz_tol)
    breaks <- c(0, indices, l)
    return(breaks)
}

#' Compute rt relative tolerance to use.
#' @description
#' Compute the elution time tolerance based on the kernel density estimation.
#' It plots the fitting function if set to TRUE. 
#' @param max.num.segments the maximum number of segments.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param number_of_samples The number of spectra in this analysis.
#' @param chr retention time of all peaks in all profiles in the study.
#' @param min.bins the minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins the maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#' @param do.plot Indicates whether plot should be drawn.
#' @return rt_tol_relative the elution time tolerance.
compute_rt_tol_relative <- function(breaks,
                                    max.num.segments,
                                    aver.bin.size,
                                    number_of_samples,
                                    chr,
                                    min.bins,
                                    max.bins,
                                    do.plot = FALSE) {
    da <- 0
    #' This conditional makes sure that length(s) is <= max.num.segments
    #' If False, length(s) =  max.num.segments, and s[i] is the largest
    #' integer no greater than the corresponding element. Otherwise
    #' length(s) =  length(breaks) - 1
    if (length(breaks) > max.num.segments) {
        s <- floor(seq(2, length(breaks), length.out = max.num.segments))
    } else {
        s <- 2:length(breaks)
    }

    #' This loop creates a vector with distances between rt peaks. Distances
    #' are stored in a triangular matrix and converted to a vector subsequently.
    #' Vector length should be < 100, otherwise, vector is 
    #' constructed extracting only 100 samples.    
    for (i in s) {
        subset_idx <- (breaks[i - 1] + 1):breaks[i]# create subset of indices
        if (length(subset_idx) <= 3 * number_of_samples) {
            rt_distances <- as.vector(dist(chr[subset_idx]))  
            if (length(rt_distances) > 100) {
                rt_distances <- sample(rt_distances, 100) 
            }
            da <- c(da, rt_distances)
        }
    }
    
    #' Calculation of kernel density estimation to estimate the rt_tol_relative
    da <- da[!is.na(da)]
    uppermost <- max(da)
    n <- min(max.bins, max(min.bins, round(length(da) / aver.bin.size)))
    des <- density(da,
        kernel = "gaussian", n = n,
        bw = uppermost / n * 2, from = 0
    )
    y <- des$y[des$x > 0]
    x <- des$x[des$x > 0]

    this.l <- lm(y[x > uppermost / 4] ~ x[x > uppermost / 4])
    exp.y <- this.l$coef[1] + this.l$coef[2] * x
    y2 <- y[1:(length(y) - 1)]
    y3 <- y[2:(length(y))]
    y2[which(y2 < y3)] <- y3[which(y2 < y3)]
    y[1:(length(y) - 1)] <- y2

    yy <- cumsum(y > 1.5 * exp.y)
    yi <- seq_along(yy)
    sel <- min(which(yy < yi)) - 1
    rt_tol_relative <- x[sel]

    if (do.plot) {
        tolerance_plot(x, y, exp.y,
            sel,
            main = "find retention time tolerance"
        )
    }
    return(rt_tol_relative)
}

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
#' @return A matrix with six columns. Every row corresponds to a peak in one of the spectrum. The columns are: m/z, elution time, spread, signal strength,
#'  spectrum label, and peak group label. The rows are ordered by the median m/z of each peak group, and with each peak group the rows are ordered
#'  by the elution time.
#' @examples
#' find.tol.time(mz, chr, lab, number_of_samples = number_of_samples, mz_tol_relative = mz_tol_relative, mz_tol_absolute = mz_tol_absolute, do.plot = FALSE)
find.tol.time <- function(mz,
                          chr,
                          lab,
                          number_of_samples,
                          mz_tol_relative = 2e-5,
                          rt_tol_relative = NA,
                          aver.bin.size = 200,
                          min.bins = 50,
                          max.bins = 100,
                          mz_tol_absolute = 0.01,
                          max.num.segments = 10000,
                          do.plot = TRUE) {
    features <- tibble::tibble(mz = mz, rt = chr, labels = lab)
    features <- dplyr::arrange_at(features, "mz")

    min_mz_tol <- compute_min_mz_tolerance(
        features$mz,
        mz_tol_relative,
        mz_tol_absolute
    )

    mz_breaks <- compute_breaks_3(features$mz, min_mz_tol)
    features$mz_group <- 0

    for (i in 2:length(mz_breaks)) {
        subset_indices <- (mz_breaks[i - 1] + 1):mz_breaks[i]
        features$mz_group[subset_indices] <- i
    }

    features <- features |> dplyr::arrange_at(c("mz_group", "rt"))

    mz_breaks <- mz_breaks[c(-1, -length(mz_breaks))]

    if (is.na(rt_tol_relative)) {
        rt_tol_relative <- compute_rt_tol_relative(
            mz_breaks,
            max.num.segments,
            aver.bin.size,
            number_of_samples,
            features$rt,
            min.bins,
            max.bins
        )
    }

    rt_diffs <- diff(features$rt)
    rt_breaks <- which(rt_diffs > rt_tol_relative)
    all.breaks <- c(0, unique(c(mz_breaks, rt_breaks)), nrow(features))
    all.breaks <- all.breaks[order(all.breaks)]

    features$grps <- 0
    for (i in 2:length(all.breaks)) {
        features$grps[(all.breaks[i - 1] + 1):all.breaks[i]] <- i
    }

    list(
        mz = features$mz,
        chr = features$rt,
        lab = features$labels,
        grps = features$grps,
        chr.tol = rt_tol_relative,
        mz.tol = mz_tol_relative
    )
}
