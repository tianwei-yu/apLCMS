#' @import dplyr foreach
NULL
#> NULL

compute_comb <- function(candi, template, this.feature, j) {
  this.comb <- dplyr::bind_rows(
    dplyr::bind_cols(candi, label = rep(template, nrow(candi))),
    dplyr::bind_cols(this.feature[, 1:2], label = rep(j, nrow(this.feature)))
  )
  this.comb <- dplyr::arrange(this.comb, this.comb[, 1])
  return(this.comb)
}

compute_sel <- function(this.comb, mz_tol_relative, rt_tol_relative) {
  l <- nrow(this.comb)
  sel <- which(this.comb[2:l, 1] - this.comb[1:(l - 1), 1] <
    mz_tol_relative * this.comb[1:(l - 1), 1] * 2 &
    abs(this.comb[2:l, 2] - this.comb[1:(l - 1), 2]) <
      rt_tol_relative & this.comb[2:l, 3] != this.comb[1:(l - 1), 3])
  return(sel)
}

compute_template_adjusted_rt <- function(this.comb, sel, j) {
  all.ftr.table <- cbind(this.comb[sel, 2], this.comb[sel + 1, 2])
  to.flip <- which(this.comb[sel, 3] == j)
  temp <- all.ftr.table[to.flip, 2]
  all.ftr.table[to.flip, 2] <- all.ftr.table[to.flip, 1]
  all.ftr.table[to.flip, 1] <- temp

  # now the first column is the template retention time.
  # the second column is the to-be-adjusted retention time

  cat(c("sample", j, "using", nrow(all.ftr.table), ", "))
  if (j %% 3 == 0) {
    cat("\n")
  }

  all.ftr.table <- all.ftr.table[order(all.ftr.table[, 2]), ]
  return(all.ftr.table)
}

compute_corrected_features <- function(this.feature, this.diff, avg_time) {
  this.feature <- this.feature[order(this.feature[, 2], this.feature[, 1]), ]
  this.corrected <- this.old <- this.feature[, 2]
  to.correct <- this.old[this.old >= min(this.diff) &
    this.old <= max(this.diff)]

  this.smooth <- ksmooth(this.diff, avg_time,
    kernel = "normal",
    bandwidth = (max(this.diff) - min(this.diff)) / 5,
    x.points = to.correct
  )

  this.corrected[between(this.old, min(this.diff), max(this.diff))] <-
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
  for (i in missing_values) {
    this.d <- abs(orig.feature[i, 2] - orig.feature[, 2])
    this.d[missing_values] <- Inf
    this.s <- which(this.d == min(this.d))[1]
    this.feature[i, 2] <- orig.feature[i, 2] + this.feature[this.s, 2] -
      orig.feature[this.s, 2]
  }
  return(this.feature)
}

add_sample_id_and_rt_cluster <- function(sample, all.ft, current_sample_id) {
    sample <- sample[order(sample[, 1], sample[, 2]), ]
    group_ids <- which(all.ft$sample_id == current_sample_id)

    sample_grouped <- cbind(all.ft$mz[group_ids], all.ft$rt[group_ids], all.ft$grps[group_ids])
    sample_grouped <- sample_grouped[order(sample_grouped[, 1], sample_grouped[, 2]), ]
    
    features <- cbind(sample, sample_id = rep(current_sample_id, nrow(sample)), cluster = sample_grouped[, 3])
    return(features)
}

#' Adjust retention time across spectra.
#'
#' This function adjusts the retention time in each LC/MS profile to achieve better between-profile agreement.
#'
#' @param extracted_features A list object. Each component is a matrix which is the output from proc.to.feature()
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the
#'  program to search for the tolerance level based on the data. This value is expressed as the
#'  percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which
#'  allows the program to search for the tolerance level based on the data.
#' @param colors The vector of colors to be used for the line plots of time adjustments. The default is NA,
#'  in which case the program uses a set of default color set.
#' @param mz_max_diff Argument passed to find.tol(). Consider only m/z diffs smaller than this value.
#'  This is only used when the mz_tol_relative is NA.
#' @param mz_tol_absolute As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable
#'  when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly
#'  influences feature matching in higher m/z range.
#' @param do.plot Indicates whether plot should be drawn.
#' @param rt_colname contains the retention time information
#' @return A list object with the exact same structure as the input object features, i.e. one matrix per profile
#'  being processed. The only difference this output object has with the input object is that the retention time
#'  column in each of the matrices is changed to new adjusted values.
#' @export
#' @examples
#' adjust.time(extracted_features, mz_max_diff = 10 * 1e-05, do.plot = FALSE)
adjust.time <- function(extracted_features,
                        mz_tol_relative = NA,
                        rt_tol_relative = NA,
                        colors = NA,
                        mz_max_diff = 1e-4,
                        mz_tol_absolute = 0.01,
                        do.plot = TRUE,
                        rt_colname = "pos") {
  number_of_samples <- nrow(summary(extracted_features))

  if (number_of_samples <= 1) {
    message("Only one sample. No need to correct for time.")
  }

  if (do.plot) {
    par(mfrow = c(2, 2))
    draw_plot(label = "Retention time \n adjustment", cex = 2)
  }

  extracted_features <- lapply(extracted_features, function(x) tibble::as_tibble(x) |> dplyr::rename(rt = pos))

  values <- concatenate_feature_tables(extracted_features, rt_colname)
  all_mz <- values$mz
  all_rt <- values$rt
  all_sample_ids <- values$sample_id

  if (is.na(mz_tol_relative)) {
    mz_tol_relative <- find.tol(
      all_mz,
      mz_max_diff = mz_max_diff,
      do.plot = do.plot
    )
  } else if (do.plot) {
    draw_plot(
      main = "m/z tolerance level given",
      label = mz_tol_relative
    )
  }

  if (!is.na(rt_tol_relative) && do.plot) {
    draw_plot(
      main = "retention time \n tolerance level given",
      label = rt_tol_relative
    )
  }

  all.ft <- find.tol.time(
    all_mz,
    all_rt,
    all_sample_ids,
    number_of_samples = number_of_samples,
    mz_tol_relative = mz_tol_relative,
    rt_tol_relative = rt_tol_relative,
    mz_tol_absolute = mz_tol_absolute,
    do.plot = do.plot
  )
  rt_tol_relative <- all.ft$rt.tol

  message("**** performing time correction ****")
  message(paste("m/z tolerance level: ", mz_tol_relative))
  message(paste("time tolerance level:", rt_tol_relative))

  for (i in 1:number_of_samples) {
    sample <- extracted_features[[i]]
    features <- add_sample_id_and_rt_cluster(
      sample,
      all.ft,
      i
    )
    extracted_features[[i]] <- features
  }

  num.ftrs <- sapply(extracted_features, nrow)#as.vector(table(all.ft$sample_id))
  template <- which(num.ftrs == max(num.ftrs))[1]
  message(paste("the template is sample", template))

  candi <- extracted_features[[template]][, 1:2]

  corrected_features <- foreach::foreach(j = 1:number_of_samples, .export = c(
    "compute_corrected_features",
    "compute_template_adjusted_rt", "compute_comb", "compute_sel"
  )) %dopar% {
    this.feature <- extracted_features[[j]]
    if (j != template) {
      this.comb <- compute_comb(candi, template, this.feature, j)

      sel <- compute_sel(this.comb, mz_tol_relative, rt_tol_relative)
      if (length(sel) < 20) {
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

    if (sum(is.na(this.feature[, 2])) > 0) {
      this.feature <- fill_missing_values(
        extracted_features[[j]],
        this.feature
      )
    }
    this.feature
  }


  if (do.plot) {
    if (is.na(colors[1])) {
      colors <- c(
        "red", "blue", "dark blue", "orange", "green", "yellow",
        "cyan", "pink", "violet", "bisque", "azure", "brown",
        "chocolate", rep("grey", number_of_samples)
      )
    }

    draw_plot(
      x = range(extracted_features[[1]][, 2]),
      y = c(-rt_tol_relative, rt_tol_relative),
      xlab = "Original Retention time",
      ylab = "Retention time deviation",
      axes = TRUE
    )

    for (i in 1:number_of_samples) {
      extracted_features[[i]] <- extracted_features[[i]][order(extracted_features[[i]][, 1], extracted_features[[i]][, 2]), ]
      points(extracted_features[[i]][, 2], corrected_features[[i]][, 2] - extracted_features[[i]][, 2],
        col = colors[i], cex = .2
      )
    }
  }

  if (exists("corrected_features")) {
    return(corrected_features)
  }
}
