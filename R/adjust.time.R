#' @import dplyr foreach
NULL
#> NULL

compute_comb <- function(template_features, features) {
  combined <- dplyr::bind_rows(
    template_features,
    dplyr::bind_cols(features |> dplyr::select(c(mz, rt)), sample_id = features$sample_id)
  )
  combined <- combined |> dplyr::arrange_at("mz")
  return(combined)
}

compute_sel <- function(combined, mz_tol_relative, rt_tol_relative) {
  l <- nrow(combined)
  sel <- which(combined$mz[2:l] - combined$mz[1:(l - 1)] <
    mz_tol_relative * combined$mz[1:(l - 1)] * 2 &
    abs(combined$rt[2:l] - combined$rt[1:(l - 1)]) <
      rt_tol_relative & combined$sample_id[2:l] != combined$sample_id[1:(l - 1)])
  return(sel)
}

compute_template_adjusted_rt <- function(combined, sel, j) {
  all_features <- cbind(combined$rt[sel], combined$rt[sel + 1])
  flip_indices <- which(combined$sample_id[sel] == j)
  temp <- all_features[flip_indices, 2]
  all_features[flip_indices, 2] <- all_features[flip_indices, 1]
  all_features[flip_indices, 1] <- temp

  # now the first column is the template retention time.
  # the second column is the to-be-adjusted retention time
  
  all_features <- all_features[order(all_features[, 2]), ]
  return(all_features)
}

compute_corrected_features <- function(features, delta_rt, avg_time) {
  features <- features[order(features$rt, features$mz), ]
  corrected <- features$rt
  original <- features$rt
  to_correct <- original[original >= min(delta_rt) &
    original <= max(delta_rt)]

  this.smooth <- ksmooth(delta_rt, avg_time,
    kernel = "normal",
    bandwidth = (max(delta_rt) - min(delta_rt)) / 5,
    x.points = to_correct
  )

  corrected[dplyr::between(original, min(delta_rt), max(delta_rt))] <-
    this.smooth$y + to_correct
  corrected[original < min(delta_rt)] <- corrected[original < min(delta_rt)] +
    mean(this.smooth$y[this.smooth$x == min(this.smooth$x)])
  corrected[original > max(delta_rt)] <- corrected[original > max(delta_rt)] +
    mean(this.smooth$y[this.smooth$x == max(this.smooth$x)])
  features$rt <- corrected
  features <- features[order(features$mz, features$rt), ]
  return(features)
}

fill_missing_values <- function(orig.feature, this.feature) {
  missing_values <- which(is.na(this.feature$rt))
  for (i in missing_values) {
    this.d <- abs(orig.feature$rt[i] - orig.feature$rt)
    this.d[missing_values] <- Inf
    this.s <- which.min(this.d)
    this.feature$rt[i] <- orig.feature$rt[i] + this.feature$rt[this.s] -
      orig.feature$rt[this.s]
  }
  return(this.feature)
}

#' @export
compute_template <- function(extracted_features) {
  num.ftrs <- sapply(extracted_features, nrow)
  template_id <- which.max(num.ftrs)
  template <- extracted_features[[template_id]]$sample_id[1]
  message(paste("the template is sample", template))

  candi <- tibble::as_tibble(extracted_features[[template_id]]) |> dplyr::select(c(mz, rt))
  template_features <- dplyr::bind_cols(candi, sample_id = rep(template, nrow(candi)))
  return(tibble::as_tibble(template_features))
}

#' @export
correct_time <- function(this.feature, template_features, mz_tol_relative, rt_tol_relative) {
    orig.features <- this.feature
    template <- unique(template_features$sample_id)[1]
    j <- unique(this.feature$sample_id)[1]

    if (j != template) {
      this.comb <- compute_comb(template_features, this.feature)
      sel <- compute_sel(this.comb, mz_tol_relative, rt_tol_relative)

      if (length(sel) < 20) {
        stop("too few, aborted")
      } else {
        all.ftr.table <- compute_template_adjusted_rt(this.comb, sel, j)
        # the to be adjusted time
        this.diff <- all.ftr.table[, 2]
        # the difference between the true time and the to-be-adjusted time
        avg_time <- all.ftr.table[, 1] - this.diff
        this.feature <- compute_corrected_features(this.feature, this.diff, avg_time)
      }
    }

    if (sum(is.na(this.feature$rt)) > 0) {
      this.feature <- fill_missing_values(
        orig.features,
        this.feature
      )
    }

  return(tibble::as_tibble(this.feature, column_name = c("mz", "rt", "sd1", "sd2", "area", "sample_id", "cluster")))
}

#' Adjust retention time across spectra.
#'
#' This function adjusts the retention time in each LC/MS profile to achieve better between-profile agreement.
#'
#' @param extracted_features A list object. Each component is a matrix which is the output from compute_clusters
#' @param mz_tol_relative The m/z tolerance level for peak alignment. This value is expressed as the
#'  percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment.
#' @param colors The vector of colors to be used for the line plots of time adjustments. The default is NA,
#'  in which case the program uses a set of default color set.
#' @param do.plot Indicates whether plot should be drawn.
#' @return A list object with the exact same structure as the input object features, i.e. one matrix per profile
#'  being processed. The only difference this output object has with the input object is that the retention time
#'  column in each of the matrices is changed to new adjusted values.
#' @export
adjust.time <- function(extracted_features,
                        mz_tol_relative,
                        rt_tol_relative,
                        colors = NA,
                        do.plot = TRUE) {
  number_of_samples <- length(extracted_features)

  if (number_of_samples <= 1) {
    message("Only one sample. No need to correct for time.")
  }

  if (do.plot) {
    par(mfrow = c(2, 2))
    draw_plot(label = "Retention time \n adjustment", cex = 2)
  }

  template_features <- compute_template(extracted_features)

  corrected_features <- foreach::foreach(features = extracted_features) %do% correct_time(
    features,
    template_features,
    mz_tol_relative,
    rt_tol_relative
  )

  if (do.plot) {
    draw_rt_correction_plot(
      colors,
      extracted_features,
      corrected_features,
      rt_tol_relative
    )
  }

  if (exists("corrected_features")) {
    return(corrected_features)
  }
}
