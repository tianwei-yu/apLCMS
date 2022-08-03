to_attach <- function(pick, number_of_samples, use = "sum") {
    strengths <- rep(0, number_of_samples)
    if (is.null(nrow(pick))) {
        # this is very strange if it can ever happen
        # maybe commas are missing? we want the same as below
        # but also why if there are no rows...
        strengths[pick[6]] <- pick[5]
        return(c(pick[1], pick[2], pick[1], pick[1], strengths))
    } else {
        for (i in seq_along(strengths)) {
            # select all areas from the same sample
            areas <- pick[pick[, 6] == i, 5]
            if (use == "sum")
                strengths[i] <- sum(areas)
            if (use == "median")
                strengths[i] <- median(areas)
            # can be NA if pick does not contain any data from a sample
        }
        # average of m/z, average of rt, min of rt, max of rt, sum/median of areas
        return(c(mean(pick[, 1]), mean(pick[, 2]), min(pick[, 1]),
                 max(pick[, 1]), strengths))
    }
}


select_rt <- function(sample, rt_tol_relative, min_occurrence, number_of_samples) {
    # Kernel Density Estimation, this time for retention time
    den <- density(sample[, 2], bw = rt_tol_relative / 1.414)
    # again select statistically significant points
    turns <- find.turn.point(den$y)
    peaks <- den$x[turns$pks]
    valleys <- den$x[turns$vlys]
    for (i in seq_along(peaks)) {
        lower_bound <- max(valleys[valleys < peaks[i]])
        upper_bound <- min(valleys[valleys > peaks[i]])
        # select data with m/z within lower and upper bound from density estimation
        selected <- which(sample[, 2] > lower_bound & sample[, 2] <= upper_bound)
        sample_grouped <- sample[selected, ]
        if (!is.null(nrow(sample_grouped))) {
            if (length(unique(sample_grouped[, 6])) >= min_occurrence) {
                # continue if data is still from at least 'min_occurrence' samples
                return(c(to_attach(sample_grouped, number_of_samples, use = "sum"),
                         to_attach(sample_grouped[, c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
                         sd(sample_grouped[, 1], na.rm = TRUE)
                       )
                )
            }
        }
    }
}


select_mz <- function(sample, mz_tol_relative, rt_tol_relative, min_occurrence, number_of_samples) {
    # Kernel Density Estimation with target standard deviation 'mz_tol_relative' multiplied by median of m/z
    des <- density(sample[, 1], bw = mz_tol_relative * median(sample[, 1]))
    # find the peaks and valleys of a density function
    turns <- find.turn.point(des$y)
    peaks <- des$x[turns$pks]
    valleys <- des$x[turns$vlys]
    for (i in seq_along(peaks)) {
        lower_bound <- max(valleys[valleys < peaks[i]])
        upper_bound <- min(valleys[valleys > peaks[i]])
        # select data with m/z within lower and upper bound from density estimation
        selected <- which(sample[, 1] > lower_bound & sample[, 1] <= upper_bound)
        sample_grouped <- sample[selected, ]
        if (!is.null(nrow(sample_grouped))) {
            # continue if data is still from at least 'min_occurrence' samples
            if (length(unique(sample_grouped[, 6])) >= min_occurrence) {
                return(select_rt(sample_grouped, rt_tol_relative, min_occurrence, number_of_samples)) 
            }
        }
    }
}


# returns a list of aligned features and original peak times
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
    
    number_of_samples <- nrow(summary(features))
    if (number_of_samples > 1) {
        values <- get_feature_values(features, rt_colname)
        mz_values <- values$mz
        rt <- values$rt
        sample_id <- values$sample_id
        
        # sort all values by m/z, if equal by rt
        ordering <- order(mz_values, rt)
        mz_values <- mz_values[ordering]
        rt <- rt[ordering]
        sample_id <- sample_id[ordering]
        
        # find relative m/z tolerance level
        if (is.na(mz_tol_relative)) {
            mz_tol_relative <- find.tol(mz_values, mz_max_diff = mz_max_diff, do.plot = do.plot)
            if (length(mz_tol_relative) == 0) {
                mz_tol_relative <- 1e-5
                warning("Automatic tolerance finding failed, 10 ppm was assigned. 
                        May need to manually assign alignment mz tolerance level.")
            }
        } else if (do.plot) {
            draw_plot(main = "alignment m/z tolerance level given", 
                      label = mz_tol_relative, cex = 1.2)
        }
        
        if (!is.na(rt_tol_relative) && do.plot) {
            draw_plot(main = "retention time \n tolerance level given", 
                      label = rt_tol_relative, cex = 1.2)
        }
        
        # find relative retention time tolerance level
        # also does some preprocessing grouping steps
        all_features <- find.tol.time(mz_values,
                                      rt,
                                      sample_id,
                                      number_of_samples = number_of_samples,
                                      mz_tol_relative = mz_tol_relative,
                                      rt_tol_relative = rt_tol_relative,
                                      mz_tol_absolute = mz_tol_absolute,
                                      do.plot = do.plot)
        rt_tol_relative <- all_features$rt.tol
        
        message("**** performing feature alignment ****")
        message(paste("m/z tolerance level: ", mz_tol_relative))
        message(paste("time tolerance level:", rt_tol_relative))
        
        # create zero vectors of length number_of_samples + 4 ?
        aligned_features <- pk.times <- rep(0, 4 + number_of_samples)
        mz_sd <- 0
        
        labels <- unique(all_features$grps)
        area <- grps <- mz_values
        
        # grouping the features based on their m/z values (assuming the tolerance level)
        sizes <- c(0, cumsum(sapply(features, nrow)))
        for (i in 1:number_of_samples) {
            sample <- features[[i]]
            # order by m/z then by rt
            sample <- sample[order(sample[, 1], sample[, 2]),]
            
            # select preprocessed features belonging to current sample
            group_ids <- which(all_features$lab == i)
            # select m/z, rt and their group ID
            sample_grouped <- cbind(all_features$mz[group_ids], all_features$rt[group_ids], all_features$grps[group_ids])
            # order them again? should be ordered already...
            sample_grouped <- sample_grouped[order(sample_grouped[, 1], sample_grouped[, 2]),]
            
            # update m/z, rt, area values with ordered ones
            mz_values[(sizes[i] + 1):sizes[i + 1]] <- sample[, 1]
            rt[(sizes[i] + 1):sizes[i + 1]] <- sample[, 2]
            area[(sizes[i] + 1):sizes[i + 1]] <- sample[, 5]
            # assign row identifier
            grps[(sizes[i] + 1):sizes[i + 1]] <- sample_grouped[, 3]
            # assign batch identifier
            sample_id[(sizes[i] + 1):sizes[i + 1]] <- i
        }
        
        # table with number of values in a group
        groups_cardinality <- table(all_features$grps)
        # count those with minimal occurrence 
        # (times 3 ? shouldn't be number of sample) !!!
        curr.row <- sum(groups_cardinality >= min_occurrence) * 3
        mz_sd <- rep(0, curr.row)
        
        sel.labels <- as.numeric(names(groups_cardinality)[groups_cardinality >= min_occurrence])
        
        # retention time alignment
        aligned_features <-
            foreach(i = seq_along(sel.labels), .combine = rbind) %do% {
                if (i %% 100 == 0)
                    gc() # call Garbage Collection for performance improvement?
                # select a group
                group_ids <- which(grps == sel.labels[i])
                if (length(group_ids) > 1) {
                    # select data from the group
                    sample <- cbind(mz_values[group_ids], rt[group_ids], rt[group_ids], rt[group_ids], area[group_ids], sample_id[group_ids])
                    # continue if data is from at least 'min_occurrence' samples
                    if (length(unique(sample[, 6])) >= min_occurrence) {
                        return(select_mz(sample, mz_tol_relative, rt_tol_relative, min_occurrence, number_of_samples))
                    }
                } else {
                    if (min_occurrence == 1) {
                        sample_grouped <- c(mz_values[group_ids], rt[group_ids], rt[group_ids], rt[group_ids], area[group_ids], sample_id[group_ids])
                        return(c(to_attach(sample_grouped, number_of_samples, use = "sum"),
                                 to_attach(sample_grouped[c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
                                 NA
                                )
                              )
                    }
                }
                return(NULL)
            }
        
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
            hist(mz_sd, xlab = "m/z SD", ylab = "Frequency",
                 main = "m/z SD distribution")
            hist(apply(pk.times[, -1:-4], 1, sd, na.rm = TRUE), 
                 xlab = "Retention time SD", ylab = "Frequency",
                 main = "Retention time SD distribution")
        }
        
        return(rec)
    } else {
        message("There is but one experiment.  What are you trying to align?")
        return(0)
    }
}
