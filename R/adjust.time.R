compute_comb <- function(candi, template, this.feature, j){
    this.comb <- dplyr::bind_rows(dplyr::bind_cols(candi, label = rep(template, nrow(candi))),
                        dplyr::bind_cols(this.feature[, 1:2], label = rep(j, nrow(this.feature))))
    this.comb <- dplyr::arrange(this.comb, this.comb[, 1])
    return(this.comb)
}

compute_sel <- function(this.comb, mz_tol_relative, rt_tol_relative){
    l <- nrow(this.comb)
    sel <- which(this.comb[2:l, 1] - this.comb[1:(l-1), 1] < 
                   mz_tol_relative * this.comb[1:(l-1), 1] * 2 & 
                   abs(this.comb[2:l, 2] - this.comb[1:(l-1), 2]) < 
                   rt_tol_relative & this.comb[2:l, 3] != this.comb[1:(l-1), 3])
    return(sel)
}

compute_template_adjusted_rt <- function(this.comb, sel, j){
    all.ftr.table <- cbind(this.comb[sel, 2], this.comb[sel+1, 2])
    to.flip <- which(this.comb[sel, 3] == j)
    temp <- all.ftr.table[to.flip, 2]
    all.ftr.table[to.flip, 2] <- all.ftr.table[to.flip, 1]
    all.ftr.table[to.flip, 1] <- temp
    
    # now the first column is the template retention time. 
    # the second column is the to-be-adjusted retention time
    
    cat(c("sample", j, "using", nrow(all.ftr.table), ", "))
    if(j %% 3 == 0) 
      cat("\n")
    
    all.ftr.table <- all.ftr.table[order(all.ftr.table[, 2]), ]
    return(all.ftr.table)
}

compute_corrected_features <- function(this.feature, this.diff, avg_time){
    this.feature <- this.feature[order(this.feature[, 2], this.feature[, 1]), ]
    this.corrected <- this.old <- this.feature[, 2]
    to.correct <- this.old[this.old >= min(this.diff) & 
                             this.old <= max(this.diff)]
    
    this.smooth <- ksmooth(this.diff, avg_time, kernel="normal", 
                           bandwidth = (max(this.diff) - min(this.diff)) / 5, 
                           x.points = to.correct)
    
    this.corrected[this.old >= min(this.diff) & this.old <= max(this.diff)] <- 
      this.smooth$y + to.correct
    this.corrected[this.old < min(this.diff)] <- this.corrected[this.old < min(this.diff)] + 
      mean(this.smooth$y[this.smooth$x == min(this.smooth$x)])
    this.corrected[this.old > max(this.diff)] <- this.corrected[this.old > max(this.diff)] + 
      mean(this.smooth$y[this.smooth$x == max(this.smooth$x)])
    this.feature[, 2] <- this.corrected
    this.feature <- this.feature[order(this.feature[, 1], this.feature[, 2]), ]
    return(this.feature)
}

fill_missing_values <- function(orig.feature, this.feature) {
    missing_values <- which(is.na(this.feature[, 2]))
    for(i in missing_values) {
        this.d <- abs(orig.feature[i, 2] - orig.feature[, 2])
        this.d[missing_values] <- Inf
        this.s <- which(this.d == min(this.d))[1]
        this.feature[i, 2] <- orig.feature[i, 2] + this.feature[this.s, 2] - 
        orig.feature[this.s, 2]
    }
    return(this.feature)
}

# features is a list project, each sub-object is a matrix as identified by prof.to.features
adjust.time <- function(features,
                        mz_tol_relative = NA,
                        rt_tol_relative = NA,
                        colors = NA,
                        mz_max_diff = 1e-4,
                        mz_tol_absolute = 0.01,
                        do.plot = TRUE,
                        rt_colname = "pos") {
    
    number_of_samples <- nrow(summary(features))
    if(number_of_samples > 1) {
        if (do.plot) {
            par(mfrow = c(2,2))
            draw_plot(label = "Retention time \n adjustment", 
                      cex = 2)
        }

        values <- get_feature_values(features, rt_colname)
        mz <- values$mz
        chr <- values$chr
        lab <- values$lab

        if(is.na(mz_tol_relative)) {
            mz_tol_relative <- find.tol(mz, mz_max_diff = mz_max_diff, do.plot = do.plot)
        } else if(do.plot) {
            draw_plot(main = "m/z tolerance level given", 
                      label = mz_tol_relative)
        }
        
        if(!is.na(rt_tol_relative) && do.plot) {
            draw_plot(main = "retention time \n tolerance level given", 
                      label = rt_tol_relative)
        }
        
        all.ft <- find.tol.time(mz,
                                chr,
                                lab,
                                number_of_samples = number_of_samples,
                                mz_tol_relative = mz_tol_relative,
                                rt_tol_relative = rt_tol_relative,
                                mz_tol_absolute = mz_tol_absolute,
                                do.plot = do.plot)
        rt_tol_relative <- all.ft$chr.tol
        
        message("**** performing time correction ****")
        message(paste("m/z tolerance level: ", mz_tol_relative))
        message(paste("time tolerance level:", rt_tol_relative))
        
        for(i in 1:number_of_samples) {
            this <- features[[i]]
            sel <- which(all.ft$lab == i)
            that <- cbind(all.ft$mz[sel], all.ft$chr[sel], all.ft$grps[sel])
            this <- this[order(this[, 1], this[, 2]), ]
            that <- that[order(that[, 1], that[, 2]), ]
            this <- cbind(this, rep(i, nrow(this)), that[, 3])
            features[[i]] <- this
        }
        
        num.ftrs <- as.vector(table(all.ft$lab))
        template <- which(num.ftrs == max(num.ftrs))[1]
        message(paste("the template is sample", template))
        
        candi <- features[[template]][, 1:2]
        
        corrected_features <- foreach(j = 1:number_of_samples,.export = c("compute_corrected_features",
        "compute_template_adjusted_rt", "compute_comb", "compute_sel")) %dopar% {
            this.feature <- features[[j]]
            if(j != template) {
                this.comb <- compute_comb(candi, template, this.feature, j)

                sel <- compute_sel(this.comb, mz_tol_relative, rt_tol_relative)
                if(length(sel) < 20) {
                    cat("too few, aborted")
                } else {
                    all.ftr.table <- compute_template_adjusted_rt(this.comb, sel, j)
                    
                    # the to be adjusted time
                    this.diff <- all.ftr.table[, 2]
                    
                    # the difference between the true time and the to-be-adjusted time
                    avg_time <- all.ftr.table[, 1] - this.diff
                    
                    this.feature <- compute_corrected_features(this.feature, this.diff, avg_time)
                }
            }

            if(sum(is.na(this.feature[, 2])) > 0) {
                this.feature <- fill_missing_values(
                    features[[j]],
                    this.feature
                )
            }
            this.feature
        }
    } else {
        message("Only one sample. No need to correct for time.")
    }

    if (do.plot) {
        if(is.na(colors[1])) 
          colors <-c ("red", "blue", "dark blue", "orange", "green", "yellow", 
                      "cyan", "pink", "violet", "bisque", "azure", "brown",
                      "chocolate", rep("grey", number_of_samples))

        draw_plot(x = range(features[[1]][, 2]), 
                  y = c(-rt_tol_relative, rt_tol_relative), 
                  xlab = "Original Retention time",
                  ylab = "Retention time deviation",
                  axes = TRUE)

        for(i in 1:number_of_samples) {
            features[[i]] <- features[[i]][order(features[[i]][, 1], features[[i]][, 2]), ]
            points(features[[i]][, 2], corrected_features[[i]][, 2] - features[[i]][, 2], 
                   col=colors[i], cex=.2)
        }
    }

    return(corrected_features)
}
