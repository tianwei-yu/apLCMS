to_attach <- function(pick, number_of_samples, use = "sum") {
    strengths <- rep(0, number_of_samples)
    if (is.null(nrow(pick))) {
        strengths[pick[6]] <- pick[5]
        return(c(pick[1], pick[2], pick[1], pick[1], strengths))
    } else {
        for (i in seq_along(strengths)) {
            if (use == "sum")
                strengths[i] <- sum(pick[pick[, 6] == i, 5])
            if (use == "median")
                strengths[i] <- median(pick[pick[, 6] == i, 5])
        }
        return(c(mean(pick[, 1]), mean(pick[, 2]), min(pick[, 1]),
                 max(pick[, 1]), strengths))
    }
}

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
#' @param rt_colname rt_colname
#' @return returns a list of aligned features and original peak times
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
    
    number_of_samples <- nrow(summary(features))
    if (number_of_samples > 1) {
        values <- get_feature_values(features, rt_colname)
        mz_values <- values$mz
        chr <- values$chr
        lab <- values$lab
        
        o <- order(mz_values, chr)
        mz_values <- mz_values[o]
        chr <- chr[o]
        lab <- lab[o]
        
        # find relative m/z tolerance level
        if (is.na(mz_tol_relative)) {
            mz_tol_relative <- find.tol(mz_values, mz_max_diff = mz_max_diff, do.plot = do.plot)
            if (length(mz_tol_relative) == 0) {
                mz_tol_relative <- 1e-5
                warning(
                    "Automatic tolerance finding failed, 10 ppm was assigned. May need to manually assign alignment mz tolerance level."
                )
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
        all.ft <- find.tol.time(mz_values,
                                chr,
                                lab,
                                number_of_samples = number_of_samples,
                                mz_tol_relative = mz_tol_relative,
                                rt_tol_relative = rt_tol_relative,
                                mz_tol_absolute = mz_tol_absolute,
                                do.plot = do.plot)
        rt_tol_relative <- all.ft$chr.tol
        
        message("**** performing feature alignment ****")
        message(paste("m/z tolerance level: ", mz_tol_relative))
        message(paste("time tolerance level:", rt_tol_relative))
        
        aligned.ftrs <- pk.times <- rep(0, 4 + number_of_samples)
        mz.sd.rec <- 0
        
        labels <- unique(all.ft$grps)
        
        area <- grps <- mz_values
        
        # grouping the features based on their m/z values (assuming the tolerance level)
        sizes <- c(0, cumsum(sapply(features, nrow)))
        for (i in 1:number_of_samples) {
            this <- features[[i]]
            sel <- which(all.ft$lab == i)
            that <- cbind(all.ft$mz[sel], all.ft$chr[sel], all.ft$grps[sel])
            this <- this[order(this[, 1], this[, 2]),]
            that <- that[order(that[, 1], that[, 2]),]
            
            mz_values[(sizes[i] + 1):sizes[i + 1]] <- this[, 1]
            chr[(sizes[i] + 1):sizes[i + 1]] <- this[, 2]
            area[(sizes[i] + 1):sizes[i + 1]] <- this[, 5]
            grps[(sizes[i] + 1):sizes[i + 1]] <- that[, 3]
            lab[(sizes[i] + 1):sizes[i + 1]] <- i
        }
        
        ttt <- table(all.ft$grps)
        curr.row <- sum(ttt >= min_occurrence) * 3
        mz.sd.rec <- rep(0, curr.row)
        
        sel.labels <- as.numeric(names(ttt)[ttt >= min_occurrence])
        
        # retention time alignment
        aligned.ftrs <-
            foreach(i = seq_along(sel.labels), .combine = rbind) %do% {
                if (i %% 100 == 0)
                    gc()
                this.return <- NULL
                sel <- which(grps == sel.labels[i])
                if (length(sel) > 1) {
                    this <- cbind(mz_values[sel], chr[sel], chr[sel], chr[sel], area[sel], lab[sel])
                    if (length(unique(this[, 6])) >= min_occurrence) {
                        this.den <- density(this[, 1], bw = mz_tol_relative * median(this[, 1]))
                        turns <- find.turn.point(this.den$y)
                        pks <- this.den$x[turns$pks]
                        vlys <- this.den$x[turns$vlys]
                        for (j in seq_along(pks)) {
                            this.lower <- max(vlys[vlys < pks[j]])
                            this.upper <- min(vlys[vlys > pks[j]])
                            this.sel <- which(this[, 1] > this.lower & this[, 1] <= this.upper)
                            that <- this[this.sel, ]
                            if (!is.null(nrow(that))) {
                                if (length(unique(that[, 6])) >= min_occurrence) {
                                    that.den <- density(that[, 2], bw = rt_tol_relative / 1.414)
                                    that.turns <- find.turn.point(that.den$y)
                                    that.pks <- that.den$x[that.turns$pks]
                                    that.vlys <- that.den$x[that.turns$vlys]
                                    for (k in seq_along(that.pks)) {
                                        that.lower <- max(that.vlys[that.vlys < that.pks[k]])
                                        that.upper <- min(that.vlys[that.vlys > that.pks[k]])
                                        thee <- that[that[, 2] > that.lower & that[, 2] <= that.upper, ]
                                        if (!is.null(nrow(thee))) {
                                            if (length(unique(thee[, 6])) >= min_occurrence) {
                                                this.return <-
                                                    c(to_attach(thee, number_of_samples, use = "sum"),
                                                      to_attach(thee[, c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
                                                      sd(thee[, 1], na.rm = TRUE)
                                                    )
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (min_occurrence == 1) {
                        thee <- c(mz_values[sel], chr[sel], chr[sel], chr[sel], area[sel], lab[sel])
                        this.return <- c(to_attach(thee, number_of_samples, use = "sum"),
                                         to_attach(thee[c(1, 2, 3, 4, 2, 6)], number_of_samples, use = "median"),
                                         NA
                        )
                    }
                }
                this.return
            }
        
        pk.times <- aligned.ftrs[, (5 + number_of_samples):(2 * (4 + number_of_samples))]
        mz.sd.rec <- aligned.ftrs[, ncol(aligned.ftrs)]
        aligned.ftrs <- aligned.ftrs[, 1:(4 + number_of_samples)]
        
        colnames(aligned.ftrs) <-
            colnames(pk.times) <- c("mz", "time", "mz.min", "mz.max", paste("exp", 1:number_of_samples))
        
        rec <- new("list")
        rec$aligned.ftrs <- aligned.ftrs
        rec$pk.times <- pk.times
        rec$mz.tol <- mz_tol_relative
        rec$chr.tol <- rt_tol_relative
        
        if (do.plot) {
            hist(mz.sd.rec, xlab = "m/z SD", ylab = "Frequency",
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
