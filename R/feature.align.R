#' @import foreach
NULL
#> NULL

to_attach <- function(peak, number_of_samples, use = "sum") {
    strengths <- rep(0, number_of_samples)
    if (is.null(nrow(peak))) {
        strengths[peak[6]] <- peak[5]
        return(c(peak[1], peak[2], peak[1], peak[1], strengths))
    } else {
        for (i in seq_along(strengths)) {
            # select all areas/RTs from the same sample
            values <- peak[peak[, 6] == i, 5]
            if (use == "sum") {
                strengths[i] <- sum(values)
            }
            if (use == "median") {
                strengths[i] <- median(values)
            }
            # can be NA if pick does not contain any data from a sample
        }
        # average of m/z, average of rt, min of m/z, max of m/z, sum/median of areas/RTs
        return(c(
            mean(peak[, 1]), mean(peak[, 2]), min(peak[, 1]),
            max(peak[, 1]), strengths
        ))
    }
}


create_output <- function(sample_grouped, number_of_samples, deviation) {
    return(c(
        to_attach(sample_grouped, number_of_samples, use = "sum"),
        to_attach(sample_grouped[, c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
        deviation #mz standard deviation
    ))

    # variable_metadata_row <- (feature_id, mean_mz, min_mz, max_mz, mean_rt, min_rt, max_rt, num_peaks, sample_presence)
    # intensity_row <- (feature_id, sample_0_intensity, sample_1_intensity, ...)
    # rt_row <- (feature_id, sample_0_rt, sample_1_rt, ...)
}


validate_contents <- function(samples, min_occurrence) {
    # validate whether data is still from at least 'min_occurrence' number of samples
    if (!is.null(nrow(samples))) {
        if (length(unique(samples[, 6])) >= min_occurrence) {
            return(TRUE)
        }
        return(FALSE)
    }
    return(FALSE)
}


find_optima <- function(data, bandwidth) {
    # Kernel Density Estimation
    den <- density(data, bw = bandwidth)
    # select statistically significant points
    turns <- find.turn.point(den$y)
    return(list(peaks = den$x[turns$pks], valleys = den$x[turns$vlys]))
}


filter_based_on_density <- function(sample, turns, index, i) {
    # select data within lower and upper bound from density estimation
    lower_bound <- max(turns$valleys[turns$valleys < turns$peaks[i]])
    upper_bound <- min(turns$valleys[turns$valleys > turns$peaks[i]])
    selected <- which(sample[, index] > lower_bound & sample[, index] <= upper_bound)
    return(sample[selected, ])
}


select_rt <- function(sample, rt_tol_relative, min_occurrence, number_of_samples) {
    # turns for rt
    turns <- find_optima(sample[, 2], bandwidth = rt_tol_relative / 1.414)
    for (i in seq_along(turns$peaks)) {
        sample_grouped <- filter_based_on_density(sample, turns, 2, i)
        if (validate_contents(sample_grouped, min_occurrence)) {
            return(create_output(sample_grouped, number_of_samples, sd(sample_grouped[, 1], na.rm = TRUE)))
        }
    }
}


select_mz <- function(sample, mz_tol_relative, rt_tol_relative, min_occurrence, number_of_samples) {
    # turns for m/z
    mz <- sample[, 1]
    turns <- find_optima(mz, bandwidth = mz_tol_relative * median(mz))
    for (i in seq_along(turns$peaks)) {
        sample_grouped <- filter_based_on_density(sample, turns, 1, i)
        if (validate_contents(sample_grouped, min_occurrence)) {
            return(select_rt(sample_grouped, rt_tol_relative, min_occurrence, number_of_samples))
        }
    }
}


create_rows <- function(features,
                        i,
                        sel.labels,
                        mz_tol_relative,
                        rt_tol_relative,
                        min_occurrence,
                        number_of_samples) {
    if (i %% 100 == 0) {
        gc()
    } # call Garbage Collection for performance improvement?

    # select a group
    group_ids <- which(features$cluster == sel.labels[i])
    if (length(group_ids) > 1) {
        # select data from the group
        # dplyr::slice(group_ids)
        
        sample <- cbind(
            features$mz[group_ids],
            features$rt[group_ids],
            features$rt[group_ids], # pseudo mz.min
            features$rt[group_ids], # pseudo mz.max
            features$area[group_ids],
            features$sample_id[group_ids]
        )
        # continue if data is from at least 'min_occurrence' samples
        if (validate_contents(sample, min_occurrence)) {
            return(select_mz(sample, mz_tol_relative, rt_tol_relative, min_occurrence, number_of_samples))
        }
    } else if (min_occurrence == 1) {
        sample_grouped <- c(
            features$mz[group_ids], features$rt[group_ids], features$rt[group_ids],
            features$rt[group_ids], features$area[group_ids], features$sample_id[group_ids]
        )
        return(create_output(sample_grouped, number_of_samples, NA))
    }
    return(NULL)
}

create_aligned_feature_table <- function(all_table,
                                         min_occurrence,
                                         number_of_samples,
                                         rt_tol_relative,
                                         mz_tol_relative) {
    aligned_features <- NULL

    # table with number of values per group
    groups_cardinality <- table(all_table$cluster)
    # count those with minimal occurrence
    # (times 3 ? shouldn't be number of samples) !!!
    sel.labels <- as.numeric(names(groups_cardinality)[groups_cardinality >= min_occurrence])

    # retention time alignment
    for (i in seq_along(sel.labels)) {
        rows <- create_rows(
            all_table,
            i,
            sel.labels,
            mz_tol_relative,
            rt_tol_relative,
            min_occurrence,
            number_of_samples
        )

        aligned_features <- rbind(aligned_features, rows)
    }
    return(aligned_features)
}

#' Align peaks from spectra into a feature table.
#'
#' Identifies which of the peaks from the profiles correspond to the same feature.
#'
#' @param features A list object. Each component is a matrix which is the output from proc.to.feature().
#' @param min_occurrence  A feature has to show up in at least this number of profiles to be included in the final result.
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the
#'  program to search for the tolerance level based on the data. This value is expressed as the
#'  percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which
#'  allows the program to search for the tolerance level based on the data.
#' @param mz_max_diff Argument passed to find.tol(). Consider only m/z diffs smaller than this value.
#'  This is only used when the mz_tol_relative is NA.
#' @param mz_tol_absolute As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable
#'  when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly
#'  influences feature matching in higher m/z range.
#' @param do.plot Indicates whether plot should be drawn.
#' @param rt_colname Name of the column containing the retention time information.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item aligned.ftrs - A matrix, with columns of m/z values, elution times, signal strengths in each spectrum.
#'   \item pk.times - A matrix, with columns of m/z, median elution time, and elution times in each spectrum.
#'   \item mz.tol - The m/z tolerance used in the alignment.
#'   \item rt.tol - The elution time tolerance in the alignment.
#' }
#' @export
#' @examples
#' feature.align(features, mz_max_diff = 10 * 1e-05, do.plot = FALSE)
feature.align <- function(features,
                          min_occurrence = 2,
                          mz_tol_relative = NA,
                          rt_tol_relative = NA,
                          mz_max_diff = 1e-4,
                          mz_tol_absolute = 0.01,
                          do.plot = TRUE,
                          rt_colname = "pos") {
    if (do.plot) {
        par(mfrow = c(3, 2))
        draw_plot(label = "Feature alignment", cex = 2)
        draw_plot()
    }

    features <- lapply(features, function(x) {
        x <- tibble::as_tibble(x)
        if ("pos" %in% colnames(x)) {
            x <- x |> dplyr::rename(rt = pos)
        }
        return(x)
    })

    number_of_samples <- length(features)
    if (number_of_samples > 1) {
        res <- compute_clusters(
            features,
            mz_tol_relative,
            mz_tol_absolute,
            mz_max_diff,
            rt_tol_relative
        )

        all_table <- dplyr::bind_rows(res$feature_tables)
        rt_tol_relative <- res$rt_tol_relative
        mz_tol_relative <- res$mz_tol_relative

        # create zero vectors of length number_of_samples + 4 ?
        aligned_features <- create_aligned_feature_table(
            all_table,
            min_occurrence,
            number_of_samples,
            rt_tol_relative,
            mz_tol_relative
        )

        # select columns: average of m/z, average of rt, min of m/z, max of m/z, median of rt per sample (the second to_attach call)
        times_start_idx <- 5 + number_of_samples
        times_end_idx <- 2 * (4 + number_of_samples)
        pk.times <- aligned_features[, times_start_idx:times_end_idx]
        # select columns: average of m/z, average of rt, min of m/z, max of m/z, sum of areas per sample (the first to_attach call)
        areas_end_idx <- 4 + number_of_samples
        aligned_features <- aligned_features[, 1:areas_end_idx]

        # rename columns on both tables, samples are called "exp_i"
        colnames(aligned_features) <-
            colnames(pk.times) <- c("mz", "time", "mz.min", "mz.max", paste("exp", 1:number_of_samples))

        # return both tables and both computed tolerances
        rec <- new("list")
        rec$aligned_features <- aligned_features
        rec$peak_times <- pk.times
        rec$mz_tol_relative <- mz_tol_relative
        rec$rt_tol_relative <- rt_tol_relative

        if (do.plot) {
            plot_rt_histograms(
                pk.times,
                aligned_features[, ncol(aligned_features)]
            )
        }

        return(rec)
    } else {
        message("There is but one experiment.  What are you trying to align?")
        return(0)
    }
}
