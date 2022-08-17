#' @import foreach
NULL
#> NULL

to_attach <- function(pick, number_of_samples, use = "sum") {
    strengths <- rep(0, number_of_samples)
    if (is.null(nrow(pick))) {
        strengths[pick[6]] <- pick[5]
        return(c(pick[1], pick[2], pick[1], pick[1], strengths))
    } else {
        for (i in seq_along(strengths)) {
            # select all areas/RTs from the same sample
            values <- pick[pick[, 6] == i, 5]
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
            mean(pick[, 1]), mean(pick[, 2]), min(pick[, 1]),
            max(pick[, 1]), strengths
        ))
    }
}


create_output <- function(sample_grouped, number_of_samples, deviation) {
    return(c(
        to_attach(sample_grouped, number_of_samples, use = "sum"),
        to_attach(sample_grouped[, c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
        deviation
    ))
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
    turns <- find_optima(sample[, 1], bandwidth = mz_tol_relative * median(sample[, 1]))
    for (i in seq_along(turns$peaks)) {
        sample_grouped <- filter_based_on_density(sample, turns, 1, i)
        if (validate_contents(sample_grouped, min_occurrence)) {
            return(select_rt(sample_grouped, rt_tol_relative, min_occurrence, number_of_samples))
        }
    }
}


create_rows <- function(i,
                        grps,
                        sel.labels,
                        mz_values,
                        rt,
                        area,
                        sample_id,
                        mz_tol_relative,
                        rt_tol_relative,
                        min_occurrence,
                        number_of_samples) {
    if (i %% 100 == 0) {
        gc()
    } # call Garbage Collection for performance improvement?
    # select a group
    group_ids <- which(grps == sel.labels[i])
    if (length(group_ids) > 1) {
        # select data from the group
        sample <- cbind(
            mz_values[group_ids], rt[group_ids], rt[group_ids],
            rt[group_ids], area[group_ids], sample_id[group_ids]
        )
        # continue if data is from at least 'min_occurrence' samples
        if (validate_contents(sample, min_occurrence)) {
            return(select_mz(sample, mz_tol_relative, rt_tol_relative, min_occurrence, number_of_samples))
        }
    } else if (min_occurrence == 1) {
        sample_grouped <- c(
            mz_values[group_ids], rt[group_ids], rt[group_ids],
            rt[group_ids], area[group_ids], sample_id[group_ids]
        )
        return(create_output(sample_grouped, number_of_samples, NA))
    }
    return(NULL)
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
#'   \item chr.tol - The elution time tolerance in the alignment.
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
        if("pos" %in% colnames(x)) {
            x <- x |> dplyr::rename(rt = pos)
        }
        return(x)
    })

    number_of_samples <- length(features)
    if (number_of_samples > 1) {
        values <- concatenate_feature_tables(features, rt_colname) |> dplyr::arrange_at(c("mz", "rt"))

        mz_values <- values$mz
        rt <- values$rt
        sample_id <- values$sample_id

        # find relative m/z tolerance level
        if (is.na(mz_tol_relative)) {
            mz_tol_relative <- find.tol(
                mz_values,
                mz_max_diff = mz_max_diff,
                do.plot = do.plot
            )
            if (length(mz_tol_relative) == 0) {
                mz_tol_relative <- 1e-5
                warning("Automatic tolerance finding failed, 10 ppm was assigned.
                        May need to manually assign alignment mz tolerance level.")
            }
        } else if (do.plot) {
            draw_plot(
                main = "alignment m/z tolerance level given",
                label = mz_tol_relative, cex = 1.2
            )
        }

        if (!is.na(rt_tol_relative) && do.plot) {
            draw_plot(
                main = "retention time \n tolerance level given",
                label = rt_tol_relative, cex = 1.2
            )
        }

        # find relative retention time tolerance level
        # also does some preprocessing grouping steps
        all_features <- find.tol.time(
            mz_values,
            rt,
            sample_id,
            number_of_samples = number_of_samples,
            mz_tol_relative = mz_tol_relative,
            rt_tol_relative = rt_tol_relative,
            mz_tol_absolute = mz_tol_absolute,
            do.plot = do.plot
        )
        rt_tol_relative <- all_features$rt.tol

        message("**** performing feature alignment ****")
        message(paste("m/z tolerance level: ", mz_tol_relative))
        message(paste("time tolerance level:", rt_tol_relative))

        # create zero vectors of length number_of_samples + 4 ?
        aligned_features <- pk.times <- NULL
        mz_sd <- 0

        labels <- unique(all_features$grps)
        area <- grps <- mz_values <- NULL

        # grouping the features based on their m/z values (assuming the tolerance level)
        sizes <- c(0, cumsum(sapply(features, nrow)))
        all_table <- NULL
        for (i in 1:number_of_samples) {
            sample <- add_sample_id_and_rt_cluster(
              features[[i]],
              all_features,
              i
            )
            all_table <- rbind(all_table, sample)
        }

        mz_values <- all_table$mz
        rt <- all_table$rt
        area <- all_table$area
        grps <- all_table$cluster
        sample_id <- all_table$sample_id

        browser()

        # table with number of values per group
        groups_cardinality <- table(grps)
        # count those with minimal occurrence
        # (times 3 ? shouldn't be number of samples) !!!
        curr.row <- sum(groups_cardinality >= min_occurrence) * 3
        mz_sd <- rep(0, curr.row)

        sel.labels <- as.numeric(names(groups_cardinality)[groups_cardinality >= min_occurrence])

        # retention time alignment
        for(i in seq_along(sel.labels)) {
            rows <- create_rows(
                i,
                grps,
                sel.labels,
                mz_values,
                rt,
                area,
                sample_id,
                mz_tol_relative,
                rt_tol_relative,
                min_occurrence,
                number_of_samples
            )

            aligned_features <- rbind(aligned_features, rows)
        }
        #aligned_features <- aligned_features[-1, ]

        # select columns: average of m/z, average of rt, min of m/z, max of m/z, median of rt per sample (the second to_attach call)
        pk.times <- aligned_features[, (5 + number_of_samples):(2 * (4 + number_of_samples))]
        mz_sd <- aligned_features[, ncol(aligned_features)]
        # select columns: average of m/z, average of rt, min of m/z, max of m/z, sum of areas per sample (the first to_attach call)
        aligned_features <- aligned_features[, 1:(4 + number_of_samples)]

        # rename columns on both tables, samples are called "exp_i"
        colnames(aligned_features) <-
            colnames(pk.times) <- c("mz", "time", "mz.min", "mz.max", paste("exp", 1:number_of_samples))

        # return both tables and both computed tolerances
        rec <- new("list")
        rec$aligned.ftrs <- aligned_features
        rec$pk.times <- pk.times
        rec$mz.tol <- mz_tol_relative
        rec$chr.tol <- rt_tol_relative

        if (do.plot) {
            hist(mz_sd,
                xlab = "m/z SD", ylab = "Frequency",
                main = "m/z SD distribution"
            )
            hist(apply(pk.times[, -1:-4], 1, sd, na.rm = TRUE),
                xlab = "Retention time SD", ylab = "Frequency",
                main = "Retention time SD distribution"
            )
        }

        return(rec)
    } else {
        message("There is but one experiment.  What are you trying to align?")
        return(0)
    }
}
