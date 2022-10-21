#' @import foreach

create_empty_tibble <- function(number_of_samples, metadata_colnames, intensity_colnames, rt_colnames) {
    features <- new("list")
    features$metadata <- tibble::as_tibble(matrix(nrow = 0, ncol = length(metadata_colnames)), .name_repair = ~ metadata_colnames)
    features$intensity <- tibble::as_tibble(matrix(nrow = 0, ncol = length(intensity_colnames)), .name_repair = ~ intensity_colnames)
    features$rt <- tibble::as_tibble(matrix(nrow = 0, ncol = length(rt_colnames)), .name_repair = ~ rt_colnames)
    return(features)
}


add_row <- function(df, data, i, column_names) {
    row <- matrix(c(i, data), nrow=1)
    colnames(row) <- column_names
    return(dplyr::bind_rows(df, tibble::as_tibble(row)))
}


create_output <- function(sample_grouped, sample_names) {
    number_of_samples <- length(sample_names)
    intensity_row <- rep(0, number_of_samples)
    rt_row <- rep(0, number_of_samples)
    sample_presence <- rep(0, number_of_samples)
    
    for (i in seq_along(intensity_row)) {
        filtered <- filter(sample_grouped, sample_id == sample_names[i])
        if (nrow(filtered) != 0) {
            sample_presence[i] <- 1
            intensity_row[i] <- sum(filtered$area)
            rt_row[i] <- median(filtered$rt)
        }
    }
    
    mz <- sample_grouped$mz
    rt <- sample_grouped$rt
    metadata_row <- c(mean(mz), min(mz), max(mz), mean(rt), min(rt), max(rt), nrow(sample_grouped), sample_presence)
    
    return(list(metadata_row = metadata_row, intensity_row = intensity_row, rt_row = rt_row))
}


validate_contents <- function(samples, min_occurrence) {
    # validate whether data is still from at least 'min_occurrence' number of samples
    if (!is.null(nrow(samples))) {
        if (length(unique(samples$sample_id)) >= min_occurrence) {
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


select_rt <- function(sample, rt_tol_relative, min_occurrence, sample_names) {
    turns <- find_optima(sample$rt, bandwidth = rt_tol_relative / 1.414)
    for (i in seq_along(turns$peaks)) {
        sample_grouped <- filter_based_on_density(sample, turns, 2, i)
        if (validate_contents(sample_grouped, min_occurrence)) {
            return(create_output(sample_grouped, sample_names))
        }
    }
}


select_mz <- function(sample, mz_tol_relative, rt_tol_relative, min_occurrence, sample_names) {
    turns <- find_optima(sample$mz, bandwidth = mz_tol_relative * median(sample$mz))
    for (i in seq_along(turns$peaks)) {
        sample_grouped <- filter_based_on_density(sample, turns, 1, i)
        if (validate_contents(sample_grouped, min_occurrence)) {
            return(select_rt(sample_grouped, rt_tol_relative, min_occurrence, sample_names))
        }
    }
}


create_rows <- function(features,
                        i,
                        sel.labels,
                        mz_tol_relative,
                        rt_tol_relative,
                        min_occurrence,
                        sample_names) {
    if (i %% 100 == 0) {
        gc()
    } # call Garbage Collection for performance improvement?
    
    sample <- dplyr::filter(features, cluster == sel.labels[i])
    if (nrow(sample) > 1) {
        if (validate_contents(sample, min_occurrence)) {
            return(select_mz(sample, mz_tol_relative, rt_tol_relative, min_occurrence, sample_names))
        }
    } else if (min_occurrence == 1) {
        return(create_output(sample_grouped, sample_names))
    }
    return(NULL)
}


create_aligned_feature_table <- function(all_table,
                                         min_occurrence,
                                         sample_names,
                                         rt_tol_relative,
                                         mz_tol_relative) {
    
    number_of_samples <- length(sample_names)
    metadata_colnames <- c("id", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", sample_names)
    intensity_colnames <- c("id", paste0(sample_names, "_intensity"))
    rt_colnames <- c("id", paste0(sample_names, "_rt"))
    
    aligned_features <- create_empty_tibble(number_of_samples, metadata_colnames, intensity_colnames, rt_colnames)
    
    # table with number of values per group
    groups_cardinality <- table(all_table$cluster)
    # count those with minimal occurrence
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
            sample_names
        )
        
        if (!is.null(rows)) {
            aligned_features$metadata <- add_row(aligned_features$metadata, rows$metadata_row, i, metadata_colnames)
            aligned_features$intensity <- add_row(aligned_features$intensity, rows$intensity_row, i, intensity_colnames)
            aligned_features$rt <- add_row(aligned_features$rt, rows$rt_row, i, rt_colnames)
        }
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
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item aligned.ftrs - A matrix, with columns of m/z values, elution times, signal strengths in each spectrum.
#'   \item pk.times - A matrix, with columns of m/z, median elution time, and elution times in each spectrum.
#'   \item mz.tol - The m/z tolerance used in the alignment.
#'   \item rt.tol - The elution time tolerance in the alignment.
#' }
#' @export
#' @examples
#' data(extracted)
#' feature.align(extracted, mz_max_diff = 10 * 1e-05, do.plot = FALSE)
feature.align <- function(features,
                          sample_names,
                          min_occurrence = 2,
                          mz_tol_relative = NA,
                          rt_tol_relative = NA,
                          mz_max_diff = 1e-4,
                          mz_tol_absolute = 0.01,
                          do.plot = TRUE) {
    if (do.plot) {
        par(mfrow = c(3, 2))
        draw_plot(label = "Feature alignment", cex = 2)
        draw_plot()
    }
    
    number_of_samples <- length(features)
    if (number_of_samples > 1) {
        res <- compute_clusters(
            sample_names,
            features,
            mz_tol_relative,
            mz_tol_absolute,
            mz_max_diff,
            rt_tol_relative
        )

        all_table <- dplyr::bind_rows(res$feature_tables)

        aligned_features <- create_aligned_feature_table(
            all_table,
            min_occurrence,
            sample_names,
            res$rt_tol_relative,
            res$mz_tol_relative
        )
        
        aligned_features$mz_tol_relative <- res$mz_tol_relative
        aligned_features$rt_tol_relative <- res$rt_tol_relative

        return(aligned_features)
    } else {
        message("There is but one experiment.  What are you trying to align?")
        return(0)
    }
}
