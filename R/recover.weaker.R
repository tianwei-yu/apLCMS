#' @import tibble dplyr
NULL
#> NULL


#' Custom way of removing duplicate rows from a specifically formatted table.
#'
#' @description
#' Rows are considered as duplicate if the 1st, 2nd and 5th column are less than 1e-10 (tolerance) apart.
#' Only a single row in this `range` is kept from a group.
#' @param new.table The table from which the duplicate rows should be removed. Needs at least 5 columns.
#' Columns 1, 2 and 5 have to be of numeric type.
#' @param tolerance Tolerance to use for numeric comparisons.
#' @return Returns the same table with duplicate rows removed.
#' @export
duplicate.row.remove <- function(features, tolerance = 1e-10) {
  new.table <- features |> dplyr::arrange_at(c("mz", "rt", "area"))
  n <- 1
  to.remove <- c()

  for (m in 2:nrow(new.table)) {
    if (abs(new.table$mz[m] - new.table$mz[n]) < tolerance &
      abs(new.table$rt[m] - new.table$rt[n]) < tolerance &
      abs(new.table$area[m] - new.table$area[n]) < tolerance) {
      to.remove <- c(to.remove, m)
    } else {
      n <- m
    }
  }

  if(length(to.remove) > 0) {
    new.table <- new.table[-to.remove, ]
  }
  new.table
}

#' Compute custom smoothed h-method derivative of function.
#' @description
#' The function adds an extrapolated element in the end and a 0 element in the front,
#' then computes the midpoints between#' neighbouring elements and then uses the `diff`
#' function to compute the changes in rt between points.
#' @param times Retention time values.
#' @return Differences between time values.
#' @export
compute_delta_rt <- function(times) {
  # add element which is 2x the last element - the second to last - basically the extrapolated next element
  all_times <- c(0, times, 2 * times[length(times)] - times[length(times) - 1])

  # basically take the mean between consecutive values as the values - somewhat smoothed
  all_times <- (all_times[1:(length(all_times) - 1)] + all_times[2:length(all_times)]) / 2

  # get the differences between the values
  all_times <- diff(all_times)
  return(all_times)
}

#' Normalize vector so that sum(vec) = 1
#' @description x / sum(x)
#' @param x Data to normalize.
#' @return Normalized data.
l2normalize <- function(x) {
  x / sum(x)
}


#' Compute the density function of mz values.
#' @description
#' The function takes the mz values and uses \link[stats]{density} to
#' compute the local density, optionally using intensity based weighting.
#' @param mz Mass values to compute the density over.
#' @param intensities Intensities of the peaks at mz values.
#' Only used if intensity_weighted == TRUE.
#' @param bandwidth Bandwidth to use to compute the kernel density.
#' @param intensity_weighted Whether to use intensity weighting or not.
#' @return \link[stats]{density} object containing the densities.
#' @export
compute_mass_density <- function(features,
                                 bandwidth,
                                 intensity_weighted,
                                 n = 512) {
  if (intensity_weighted) {
    weights <- l2normalize(features$intensities)
  } else {
    weights <- NULL
  }
  mass_density <- density(
    features$mz,
    weights = weights,
    bw = bandwidth,
    n = n
  )
  return(mass_density)
}

#' Compute custom chromatographic tolerance.
#' @description
#' Compute chromatographic tolerance for each feature. If `use_observed_range == TRUE`,
#' the whole range of retention times for all peaks is used to compute the tolerance,
#' otherwise `rt_range` is used for each feature.
#' @param use_observed_range bool Whether to use the observed chromatographic range for computation or not.
#' @param peak_rts data.frame Retention time cross table with all peak rts.
#' @param rt_range float Default chromatographic tolerance to use.
#' @param aligned_features data.frame Aligned feature table.
#' @return vector Custom chromatographic tolerances to use for each feature.
#' @export
get_custom_rt_tol <- function(use_observed_range,
                               peak_rts,
                               rt_range,
                               aligned_features) {
  custom_rt_tol <- rep(rt_range, nrow(aligned_features))

  if (use_observed_range) {
    # check observed rt range across ALL SAMPLES
    all_peak_rts <- peak_rts[, 5:ncol(peak_rts)]
    observed.rt.range <- (apply(all_peak_rts, 1, max) - apply(all_peak_rts, 1, min)) / 2
    sufficient_rts <- apply(!is.na(all_peak_rts), 1, sum) >= 5
    selection <- which(sufficient_rts & custom_rt_tol > observed.rt.range)
    custom_rt_tol[selection] <- observed.rt.range[selection]
  }

  return(custom_rt_tol)
}

#' Compute target times for regions of interest for recovery.
#' @description
#' Compute the individual target times for the features to be recovered in the sample.
#' Spline interpolation using \link[splines]{interpSpline} is used to map from adjusted times
#' back into the original times. The function requires `x` to be distinct, hence the filtering
#' to only include rt values that occurr only once in both lists.
#' @param aligned_rts vector Aligned retention time values.
#' @param original_sample_rts vector Original feature retention time values before correction.
#' @param adjusted_sample_rts vector Feature retention time values after time correction.
#' @export
compute_target_times <- function(aligned_rts,
                                 original_sample_rts,
                                 adjusted_sample_rts,
                                 min_common_times = 4) {
  to.use <- get_times_to_use(original_sample_rts, adjusted_sample_rts)
  adjusted_subset <- adjusted_sample_rts[to.use]

  sel_non_na <- which(!is.na(aligned_rts))
  if (length(adjusted_subset) >= min_common_times) {
    original_subset <- original_sample_rts[to.use]
    sp <- splines::interpSpline(
      original_subset ~ adjusted_subset,
      na.action = na.omit
    )
    aligned_rts[sel_non_na] <- predict(sp, aligned_rts[sel_non_na])$y
  }
}

#' Get boolean mask for values that occur only once.
#' @description
#' Uses the \link[base]{table} function to compute the occurrences and then
#' checks which values only occur a single time.
#' @param values vector Values for which to compute the mask.
#' @return vector Boolean vector which is the mask of values occuring only once.
get_single_occurrence_mask <- function(values) {
  ttt <- table(values)
  mask <- values %in% as.numeric(names(ttt)[ttt == 1])
  return(mask)
}

#' Get retention time values to use
#' @description
#' Obtain retention time values which occur only once in both the original and the adjusted times.
#' This is a custom version of the unique or intersection function with rounding etc.
#' @param original_sample_rts vector Original feature retention time values before correction.
#' @param adjusted_sample_rts vector Feature retention time values after time correction.
#' @param cap int Maximum number of time values to return.
#' @return Indices of retention time values to use.
#' @export
get_times_to_use <- function(original_sample_rts, adjusted_sample_rts, cap = 2000) {
  to.use <- which(
    get_single_occurrence_mask(adjusted_sample_rts) &
      get_single_occurrence_mask(original_sample_rts)
  )

  if (length(to.use) > cap) {
    to.use <- sample(to.use, cap, replace = FALSE)
  }
  return(to.use)
}

#' Predict the indices for the valley points with low mass density.
#' @description
#' The density of mz values in the feature table is computed based on the tolerance.
#' The valleys or breaks of clusters in mz values are detected and a function
#' is approximated to predict the indices for the mass values which are the closest to those
#' valley points.
#' @param features data.frame Data table with features for which to predict the indices.
#' @param mz_orig_tol float Mz tolerance to use as KDE bandwidth parameter.
#' @return vector Predicted indices for valley points.
#' @export
predict_mz_break_indices <- function(features, mz_orig_tol) {
  mz_density <- compute_mass_density(
    features,
    TRUE,
    bandwidth = 0.5 * mz_orig_tol * max(features$mz),
    n = 2^min(15, floor(log2(length(features$mz))) - 2)
  )

  turnpoints <- find.turn.point(mz_density$y)
  mz_valleys <- mz_density$x[turnpoints$vlys]

  indices <- seq_along(features$mz)
  predictions <- approx(
    x = features$mz,
    y = indices,
    xout = mz_valleys,
    rule = 2,
    ties = "ordered"
  )$y

  predicted_indices_at_vlys <- unique(round(predictions))[-1]
  breaks <- c(0, predicted_indices_at_vlys)
  return(breaks)
}

#' Compute range of valley indices which are in mz_tol range around aligned_feature_mass.
#' @description
#'
#' @param aligned_feature_mass float Mz value of the aligned feature.
#' @param mz vector mz values of the features.
#' @param vlys_indices vector Indices of the valley points of mz clusters.
#' @param mz_tol float Tolerance to use to check if values are close.
#' @return pair Index range (start, end).
#' @export
get_mzrange_bound_indices <- function(aligned_feature_mass,
                                      mz,
                                      vlys_indices,
                                      mz_tol) {
  if (aligned_feature_mass <= mz[vlys_indices[2]]) {
    all_indices <- c(1, 2)
  } else {
    # get all indices where mz is close to aligned_feature mass,
    # the first one that is larger and the last one that is smaller.
    upper_bound_idx <- min(which(mz[vlys_indices] > aligned_feature_mass))
    lower_bound_idx <- max(which(mz[vlys_indices] < aligned_feature_mass))
    valley_indices_within_tol <- which(abs(mz[vlys_indices] - aligned_feature_mass) < mz_tol)
    all_indices <- c(
      valley_indices_within_tol,
      upper_bound_idx,
      lower_bound_idx
    ) + 1
  }
  return(list(start = min(all_indices), end = max(all_indices)))
}

#' Get indices where rt in `features` is within `rt_tol` of `target_time`.
#' @param target_time float Target retention time region.
#' @param features tibble Feature table including `labels` column.
#' @param rt_tol float Retention time tolerance.
#' @return vector Indices which are within `rt_tol` from `target_time` or
#' 1 if `target_time` is NA.
#' @export
get_rt_region_indices <- function(target_time, features, rt_tol) {
  if (!is.null(target_time) && !is.na(target_time)) {
    selection <- which(abs(features$labels - target_time) < rt_tol)
  } else {
    selection <- 1:nrow(features)
  }
  return(selection)
}

#' Get peaks and valleys of smoothed rt values in range.
#'
#' @param features tibble Data table with `labels` and `intensities` columns.
#' @param times vector Raw retention time data from raw data file.
#' @param bw float Bandwidth to use for kernel smoothing.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item pks - vector - The data points at which the density peaks.
#'   \item vlys - vector - The points in the data where the density is low
#'                         (forming a valley in the function).
get_features_in_rt_range <- function(features, times, bw) {
  time_curve <- times[between(times, min(features$labels), max(features$labels))]

  this.curve <- cbind(time_curve, time_curve * 0)
  this.curve[this.curve[, 1] %in% features$labels, 2] <- features$intensities

  this.smooth <- ksmooth(
    this.curve[, 1],
    this.curve[, 2],
    kernel = "normal",
    bandwidth = bw
  )

  return(compute_peaks_and_valleys(this.smooth))
}

#' Count the number of peaks in all valleys
#' @description
#' For each peak in ROI, count the peaks between the valley points.
#' @param roi list Named list with vectors `pks` and `vlys`.
#' @param times vector Retention time values
#' @return vector Numbers of peaks within each region defined by a peak and the two valley points.
count_peaks <- function(roi, times) {
  num_peaks <- rep(0, length(roi$pks))

  for (m in seq_along(roi$pks)) {
    boundaries <- compute_boundaries(roi$vlys, roi$pks[m])
    num_peaks[m] <- sum(between(times, boundaries$lower, boundaries$upper))
  }
  return(num_peaks)
}

#' Compute peaks and valleys which have at least `recover_min_count` peaks.
#'
#' @param features tibble Features with `mz`, `labels` and `intensities`.
#' @param times vector Retention time values from the raw data file.
#' @param bandwidth float Bandwidth to use in smoothing.
#' @param target_rt float Retention time at which to recover the intensity.
#' @param recover_min_count int Minimum number of peaks required in the area to recover the signal.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item pks - vector - The data points at which the density peaks with at least `recover_min_count` peaks between the valley points.
#'   \item vlys - vector - The points in the data where the density is low
#'                         (forming a valley in the function).
compute_pks_vlys_rt <- function(features, times, bandwidth, target_rt, recover_min_count) {
  roi <- get_features_in_rt_range(
    features,
    times,
    bandwidth
  )

  pks <- roi$pks
  vlys <- roi$vlys

  num_peaks <- count_peaks(roi, features$labels)

  if (!is.null (target_rt) && !is.na(target_rt)) {
    pks.d <- abs(pks - target_rt) # distance from the target peak location
    pks.d[num_peaks == 0] <- Inf
    pks <- pks[which.min(pks.d)]
  } else {
    pks <- pks[num_peaks > recover_min_count]
  }

  return(list(pks = pks, vlys = vlys))
}

#' Compute interpolated retention time, its standard deviation, and intensity values,.
#'
#' @param features tibble Features with `labels` and `intensities` columns.
#' @param aver_diff float Average retention time difference.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item intensity - float - Interpolated intensity value.
#'   \item label - float - Interpolated retention time value.
#'   \item sigma - float - Standard deviation of retention times
#' @export
compute_mu_sc_std <- function(features, aver_diff) {
  x <- features$labels
  y <- features$intensities

  sum_y <- sum(y)
  miu <- sum(x * y) / sum_y # weighted retention time values
  sigma <- sqrt(sum(y * (x - miu)^2) / sum_y)
  if (sigma == 0) {
    sc <- sum_y * aver_diff
    miu <- miu
  } else {
    fitted <- dnorm(x, mean = miu, sd = sigma)
    selection <- y > 0 & fitted / dnorm(miu, mean = miu, sd = sigma) > 1e-2
    sc <- exp(sum(fitted[selection]^2 * log(y[selection] / fitted[selection]) / sum(fitted[selection]^2)))
  }

  return(list(intensity = sc, label = miu, sigma = sigma))
}

#' Compute the rectangle around recovered features given that enough peaks are present.
#' @description
#'
#' @param mz Mz value of the feature.
#' @param peak Peak around which to detect the new feature.
#' @param valleys Valley points to compute the boundary region.
#' @param features tibble Tibble with `labels` and `intensities` column.
#' @param aver_diff float Average retention time difference.
#' @param times vector Raw retention time values from raw data file.
#' @param delta_rt vector Differences between consecutive retention time values (diff(times)).
#' @importFrom dplyr between
#' @return list Triplet of mz, label and intensity for the feature.
compute_curr_rec_with_enough_peaks <- function(mz,
                                               peak,
                                               valleys,
                                               features,
                                               aver_diff,
                                               times,
                                               delta_rt) {
  # same filtering of peaks as in compute_pks_vlyws and as above
  boundaries <- compute_boundaries(valleys, peak)

  subset <- features |>
    dplyr::filter(between(labels, boundaries$lower, boundaries$upper))

  if (nrow(subset) == 1) {
    intensity <- subset$intensities * aver_diff
    label <- subset$labels
  } else if (nrow(subset) >= 10) {
    res <- compute_mu_sc_std(subset, aver_diff)
    intensity <- res$intensity
    label <- res$label
  } else {
    intensity <- interpol.area(
      subset$labels,
      subset$intensities,
      times,
      delta_rt
    )
    label <- median(subset$labels)
  }

  return(c(mz, label, intensity))
}

#' Compute bounds of area using given peak and mass valley points.
#' @description
#' The lower bound is the mass of the valley the closest but smaller than peak
#' and the upper bound is the mass of the valley the closest but higher than
#' the peak.
#' @param valley_points vector values of valley points defining clusters.
#' @param peak double Value of the peak for which to get the valley bounds.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item lower - double - The value of the lower bound valley point.
#'   \item upper - double - The value of the upper bound valley point.
#' }
#' @export
compute_boundaries <- function(valley_points, peak) {
  lower <- max(valley_points[valley_points < peak])
  upper <- min(valley_points[valley_points > peak])
  return(list(lower = lower, upper = upper))
}

#' Compute peaks and valleys of density function.
#' @description
#' Given a density function, the turn points are computed and
#' the peaks and valleys in the original data (points with highest
#' and lowest density) are returned.
#' @param density stats::density Density object for which to compute peaks
#' and valleys.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item pks - vector - The data points at which the density peaks.
#'   \item vlys - vector - The points in the data where the density is low
#'                         (forming a valley in the function).
#' }
compute_peaks_and_valleys <- function(dens) {
  turns <- find.turn.point(dens$y)
  pks <- dens$x[turns$pks] # mz values with highest density
  vlys <- dens$x[turns$vlys]
  vlys <- c(-Inf, vlys, Inf) # masses with lowest densities values -> valley
  return(list(pks = pks, vlys = vlys))
}

#' Compute rectangle around feature with `aligned_feature_mz` and `target_rt` for recovery.
#'
#' @param data_table tibble Feature table with `mz`, `labels` and `intensities` column.
#' @param aligned_feature_mz float Mz value of feature in aligned feature table.
#' @param breaks vector Integer boundaries of clusters in mz values.
#' @param custom_mz_tol float Custom mz tolerance for the feature.
#' @param orig_mz_tol float Flat original mz tolerance to use.
#' @param use_intensity_weighting bool Whether to use intensity weighting.
#' @param recover_min_count int Minimum number of peaks required in the area to recover the signal.
#' @param target_rt float Target retention time value.
#' @param custom_rt_tol float Custom chromatographic tolerance to use.
#' @param times vector Raw retention time values from raw data file.
#' @param delta_rt vector Differences between consecutive retention time values (diff(times)).
#' @param aver_diff float Average retention time difference.
#' @param bandwidth float Bandwidth to use in smoothing.
#' @param min.bw float Minimum bandwidth to use.
#' @param max.bw float Maximum bandwidth to use.
#' @return tibble Tibble with `mz`, `labels` and `intensities` columns.
compute_rectangle <- function(data_table,
                              aligned_feature_mz,
                              breaks,
                              custom_mz_tol,
                              orig_mz_tol,
                              use_intensity_weighting,
                              recover_min_count,
                              target_rt,
                              custom_rt_tol,
                              times,
                              delta_rt,
                              aver_diff,
                              bandwidth,
                              min.bw,
                              max.bw) {
  bounds <- get_mzrange_bound_indices(
    aligned_feature_mz,
    data_table$mz,
    breaks,
    orig_mz_tol
  )

  features <- dplyr::slice(
    data_table,
    (breaks[bounds$start] + 1):breaks[bounds$end]
  ) |> dplyr::arrange_at("labels")

  mass.den <- compute_mass_density(
    features,
    bandwidth = 0.5 * orig_mz_tol * aligned_feature_mz,
    intensity_weighted = use_intensity_weighting
  )

  # find peaks in mz range in raw data
  mass_range <- compute_peaks_and_valleys(mass.den)
  pks_in_tol <- abs(mass_range$pks - aligned_feature_mz) < custom_mz_tol / 1.5
  mass_range$pks <- mass_range$pks[pks_in_tol]

  this.rec <- tibble::tibble(mz = Inf, labels = Inf, intensities = Inf)
  for (peak in mass_range$pks) {
    # get mass values of valleys the closest to the peak
    mass <- compute_boundaries(mass_range$vlys, peak)

    that <- features |> dplyr::filter(mz > mass$lower & mz <= mass$upper)

    # get values in RT region of interest?
    if (nrow(that) > recover_min_count) {
      that.prof <- combine.seq.3(that)
      that.mass <- sum(that.prof$mz * that.prof$intensities) / sum(that.prof$intensities)
      curr.rec <- c(that.mass, NA, NA)

      if (nrow(that.prof) < 10) {
        thee.sel <- get_rt_region_indices(
          target_rt,
          that.prof,
          custom_rt_tol
        )

        if (length(thee.sel) > recover_min_count) {
          if (length(thee.sel) > 1) {
            curr.rec[3] <- interpol.area(
              that.prof$labels[thee.sel],
              that.prof$intensities[thee.sel],
              times,
              delta_rt
            )
          } else {
            curr.rec[3] <- that.prof$intensities[thee.sel] * aver_diff
          }
          curr.rec[2] <- median(that.prof$labels[thee.sel])
          this.rec <- tibble::add_row(
            this.rec,
            tibble::tibble_row(
              mz = curr.rec[1],
              labels = curr.rec[2],
              intensities = curr.rec[3]
            )
          )
        }
      } else {
        labels_intensities <- dplyr::select(
          that.prof,
          c("labels", "intensities")
        ) |> dplyr::arrange_at("labels")
        bw <- min(max(bandwidth * (span(labels_intensities$labels)), min.bw), max.bw)

        all <- compute_pks_vlys_rt(
          labels_intensities,
          times,
          bw,
          target_rt,
          recover_min_count
        )

        for (peak in all$pks) {
          curr.rec <- compute_curr_rec_with_enough_peaks(
            that.mass,
            peak,
            all$vlys,
            labels_intensities,
            aver_diff,
            times,
            delta_rt
          )

          this.rec <- tibble::add_row(
            this.rec,
            tibble::tibble_row(
              mz = curr.rec[1],
              labels = curr.rec[2],
              intensities = curr.rec[3]
            )
          )
        }
      }
    }
  }
  return(this.rec)
}

#' Refine the selection based on mz and rt differences.
#' @param target_rt float Target retention time value.
#' @param rectangle tibble Features with columns `labels` and `mz`.
#' @param aligned_mz float Mz value in the aligned feature table of the
#' feature to be recovered.
#' @param rt_tol float Retention time tolerance.
#' @param mz_tol float Mz tolerance to use.
#' @return int Index of value in rectable closest to `target_rt` and `aligned_mz`.
refine_selection <- function(target_rt, rectangle, aligned_mz, rt_tol, mz_tol) {
  if (!is.na(target_rt)) {
    rt_term <- (rectangle$labels - target_rt)^2 / rt_tol^2
    mz_term <- (rectangle$mz - aligned_mz)^2 / mz_tol^2
    this.d <- rt_term + mz_term
  } else {
    this.d <- abs(rectangle$mz - aligned_mz)
  }
  this.sel <- which.min(this.d)
  return(this.sel)
}

#' Recover weak signals in some profiles that is not identified as a peak, but corresponds to identified peaks in other spectra.
#'
#' @description
#' Given the aligned feature table, some features are identified in a subgroup of spectra. This doesn't mean they don't exist in the other spectra.
#' The signal could be too low to pass the run filter. Thus after obtaining the aligned feature table, this function re-analyzes each spectrum to
#' try and fill in the holes in the aligned feature table.
#' @param filename the cdf file name from which weaker signal is to be recovered.
#' @param loc the location of the filename in the vector of filenames.
#' @param aligned.ftrs matrix, with columns of m/z values, elution times, signal strengths in each spectrum.
#' @param pk.times matrix, with columns of m/z, median elution time, and elution times in each spectrum.
#' @param align.mz.tol the m/z tolerance used in the alignment.
#' @param align.rt.tol the elution time tolerance in the alignment.
#' @param extracted_features The matrix which is the output from proc.to.feature().
#' @param adjusted_features The matrix which is the output from proc.to.feature().
#' The retention time in this object have been adjusted by the function adjust.time().
#' @param mz.range The m/z around the feature m/z to search for observations.
#' The default value is NA, in which case 1.5 times the m/z tolerance in the aligned object will be used.
#' @param rt.range The retention time around the feature retention time to search for observations.
#' The default value is NA, in which case 0.5 times the retention time tolerance in the aligned object will be used.
#' @param use.observed.range If the value is TRUE, the actual range of the observed locations
#' of the feature in all the spectra will be used.
#' @param orig.tol The mz.tol parameter provided to the proc.cdf() function. This helps retrieve the intermediate file.
#' @param min.bw The minimum bandwidth to use in the kernel smoother.
#' @param max.bw The maximum bandwidth to use in the kernel smoother.
#' @param bandwidth A value between zero and one. Multiplying this value to the length of the signal along the
#' time axis helps determine the bandwidth in the kernel smoother used for peak identification.
#' @param recover.min.count minimum number of raw data points to support a recovery.
#' @param intensity.weighted Whether to use intensity to weight mass density estimation.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item aligned.ftrs - A matrix, with columns of m/z values, elution times, and signal strengths in each spectrum.
#'   \item pk.times - A matrix, with columns of m/z, median elution time, and elution times in each spectrum.
#'   \item mz.tol - The m/z tolerance in the aligned object.
#'   \item rt.tol - The elution time tolerance in the aligned object.
#' }
#' @export
#' @examples
#' recover.weaker(filename, loc, aligned.ftrs, pk.times, align.mz.tol, align.rt.tol, this.f1, this.f2)
recover.weaker <- function(filename,
                           sample_name,
                           aligned.ftrs,
                           pk.times,
                           align.mz.tol,
                           align.rt.tol,
                           extracted_features,
                           adjusted_features,
                           mz.range = NA,
                           rt.range = NA,
                           use.observed.range = TRUE,
                           orig.tol = 1e-5,
                           min.bw = NA,
                           max.bw = NA,
                           bandwidth = .5,
                           recover.min.count = 3,
                           intensity.weighted = FALSE) {

  # load raw data
  this.raw <- load_file(filename)
  times <- this.raw$times
  data_table <- tibble::tibble(
    mz = this.raw$masses,
    labels = this.raw$labels,
    intensities = this.raw$intensi
  ) |> dplyr::arrange_at("mz")
  rm(this.raw)

  # Initialize parameters with default values
  if (is.na(mz.range)) mz.range <- 1.5 * align.mz.tol
  if (is.na(rt.range)) rt.range <- align.rt.tol / 2
  if (is.na(min.bw)) min.bw <- span(times) / 60
  if (is.na(max.bw)) max.bw <- span(times) / 15
  if (min.bw >= max.bw) min.bw <- max.bw / 4


  times <- sort(unique(times))
  aver.diff <- mean(diff(times))
  vec_delta_rt <- compute_delta_rt(times)

  sample_intensities <- aligned.ftrs[, sample_name]
  sample_times <- pk.times[, sample_name]

  custom.mz.tol <- mz.range * aligned.ftrs$mz
  custom.rt.tol <- get_custom_rt_tol(
    use.observed.range,
    pk.times,
    rt.range,
    aligned.ftrs
  )

  # # rounding is used to create a histogram of retention time values
  target_times <- compute_target_times(
    aligned.ftrs[, "rt"],
    round(extracted_features$pos, 5),
    round(adjusted_features$rt, 5)
  )

  # IMPORTANT: THIS CODE SECTION COULD BE USED TO REPLACE COMPUTE_TARGET_TIMES FOR THE TEST CASES AND
  # IS A MASSIVE SIMPLIFICATION.
  # sp <- splines::interpSpline(
  #   unique(extracted_features$pos) ~ unique(adjusted_features$rt),
  #   na.action = na.omit
  # )
  # target_times <- predict(sp, aligned.ftrs[, "rt"])$y


  breaks <- predict_mz_break_indices(data_table, orig.tol)

  this.mz <- rep(NA, length(sample_intensities))
  max_mz <- max(data_table$mz)

  # THIS CONSTRUCT TO EXTRACT MISSING FEATURES COULD BE USED TO POSSIBLY SPEED UP
  # THE COMPUTATION AS THE LOOP WILL ONLY GO OVER THE ROWS AND THE ADDITIONAL VARIABLES
  # CAN BE ADDED USING THE MUTATE FUNCTION.
  # t_aligned <- tibble::tibble(aligned.ftrs)
  # missing_features <- dplyr::filter(t_aligned, !!rlang::sym(sample_name) == 0 & mz < max_mz)

  # if(nrow(missing_features) > 0) {
  #   browser()
  # }

  for (i in seq_along(sample_intensities))
  {
    if (sample_intensities[i] == 0 && aligned.ftrs[i, "mz"] < max_mz) {
      this.rec <- compute_rectangle(
        data_table,
        aligned.ftrs[i, "mz"],
        breaks,
        custom.mz.tol[i],
        orig.tol,
        intensity.weighted,
        recover.min.count,
        target_times[i],
        custom.rt.tol[i] * 2,
        times,
        vec_delta_rt,
        aver.diff,
        bandwidth,
        min.bw,
        max.bw
      )

      this.sel <- get_rt_region_indices(
        target_times[i],
        this.rec,
        custom.rt.tol[i]
      )
      this.sel <- this.sel[this.sel != 1]

      if (length(this.sel) > 0) {
        if (length(this.sel) > 1) {
          this.sel <- refine_selection(
            target_times[i],
            this.rec,
            aligned.ftrs[i, 1],
            custom.rt.tol[i],
            custom.mz.tol[i]
          )
        }

        this.pos.diff <- which.min(abs(extracted_features$pos - this.rec$labels[this.sel]))
        extracted_features <- extracted_features |> tibble::add_row(
          mz = this.rec$mz[this.sel],
          pos = this.rec$labels[this.sel],
          area = this.rec$intensities[this.sel]
        )
        
        this.time.adjust <- (-extracted_features$pos[this.pos.diff] + adjusted_features$rt[this.pos.diff])

        adjusted_features <- adjusted_features |> tibble::add_row(
          mz = this.rec$mz[this.sel],
          rt = this.rec$labels[this.sel] + this.time.adjust,
          area = this.rec$intensities[this.sel],
          sample_id = grep(sample_name, colnames(aligned.ftrs)) - 4 # offset for other columns `mz`, `rt` etc
        )

        sample_intensities[i] <- this.rec$intensities[this.sel]
        sample_times[i] <- this.rec$labels[this.sel] + this.time.adjust
        this.mz[i] <- this.rec$mz[this.sel]
      }
    }
  }
  to.return <- new("list")
  to.return$this.mz <- this.mz
  to.return$this.ftrs <- sample_intensities
  to.return$this.times <- sample_times
  to.return$this.f1 <- duplicate.row.remove(extracted_features |> dplyr::rename(rt = pos)) |> dplyr::rename(pos = rt)
  to.return$this.f2 <- duplicate.row.remove(adjusted_features)

  return(to.return)
}
