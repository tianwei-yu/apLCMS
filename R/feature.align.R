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


create_output <- function(sample_grouped, number_of_samples, deviation) {
    return(c(to_attach(sample_grouped, number_of_samples, use = "sum"),
             to_attach(sample_grouped[, c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
             deviation
    )
    )
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
        
        # table with number of values per group
        groups_cardinality <- table(all_features$grps)
        # count those with minimal occurrence 
        # (times 3 ? shouldn't be number of samples) !!!
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
                    sample <- cbind(mz_values[group_ids], rt[group_ids], rt[group_ids], 
                                    rt[group_ids], area[group_ids], sample_id[group_ids])
                    # continue if data is from at least 'min_occurrence' samples
                    if (validate_contents(sample, min_occurrence)) {
                        return(select_mz(sample, mz_tol_relative, rt_tol_relative, min_occurrence, number_of_samples))
                    }
                } else if (min_occurrence == 1) {
                    sample_grouped <- c(mz_values[group_ids], rt[group_ids], rt[group_ids], 
                                        rt[group_ids], area[group_ids], sample_id[group_ids])
                    return(create_output(sample_grouped, number_of_samples, NA))
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
