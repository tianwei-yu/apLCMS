#' @import parallel doParallel
NULL
#> NULL

.merge_peaks <- function(aligned, known_table, match_tol_ppm) {
  if (is.na(match_tol_ppm)) {
    match_tol_ppm <- aligned$mz_tolerance * 1e+06
  }

  features <- tibble::as_tibble(aligned$int_crosstab)
  known_mz <- known_table[, 6]
  known_rt <- known_table[, 11]

  mass_d2 <- mass.match(features$mz, known_mz, match_tol_ppm)
  mass_matched_pos <- which(mass_d2 > 0)

  known_assigned <- rep(0, nrow(known_table))
  new_assigned <- rep(0, nrow(features))
  pairing <- matrix(0, ncol = 2, nrow = 0)
  colnames(pairing) <- c("new", "known")

  for (i in mass_matched_pos) {
    if (new_assigned[i] == 0) {
      # find all potentially related known/newly found peaks
      prev_sel_new <- i
      threshold <- features$mz[i] * match_tol_ppm / 1e+06

      sel_known <- which(abs(known_mz - features$mz[i]) < threshold)
      sel_new <- NULL
      for (m in seq_along(sel_known)) {
        distance <- abs(features$mz - known_mz[sel_known[m]])
        sel_new <- c(sel_new, which(distance < threshold))
      }
      sel_known <- unique(sel_known)
      sel_new <- unique(sel_new)

      while (length(sel_new) > length(prev_sel_new)) {
        prev_sel_new <- sel_new

        sel_known <- NULL
        for (m in seq_along(sel_new)) {
          distance <- abs(known_mz - features$mz[sel_new[m]])
          sel_known <- c(sel_known, which(distance < threshold))
        }

        sel_new <- NULL
        for (m in seq_along(sel_known)) {
          distance <- abs(features$mz - known_mz[sel_known[m]])
          sel_new <- c(sel_new, which(distance < threshold))
        }

        sel_known <- unique(sel_known)
        sel_new <- unique(sel_new)
      }

      time_matched <- mass_matched <-
        matrix(data = 0, nrow = length(sel_known), ncol = length(sel_new))

      for (k in seq_along(sel_known)) {
        time_matched[k, ] <- abs(features$rt[sel_new] - known_rt[sel_known[k]])
        mass_matched[k, ] <- abs(features$mz[sel_new] - known_mz[sel_known[k]])
      }
      mass_matched <- mass_matched/median(known_mz[sel_known])
      time_matched[mass_matched <= match_tol_ppm * 1e-06] <- 1e+10

      time_matched[is.na(time_matched)] <- aligned$rt_tolerance / 2
      both_matched <- find.match(time_matched, aligned$rt_tolerance / 2)

      for (m in seq_along(sel_new)) {
        k <- which(both_matched[, m] == 1)

        if (length(k) == 1 && known_assigned[sel_known[k]] == 0) {
          new_assigned[sel_new[m]] <- 1
          known_assigned[sel_known[k]] <- 1
          pairing <- rbind(pairing, c(sel_new[m], sel_known[k]))
        }
      }
    }
  }

  pairing
}

#' @export
augment_with_known_features <- function(aligned, known_table, match_tol_ppm) {
  pairing <- .merge_peaks(aligned, known_table, match_tol_ppm)

  features <- aligned$int_crosstab
  n_entries <- nrow(known_table) - nrow(pairing)
  to_add_ftrs <- matrix(0, ncol = ncol(features), nrow = n_entries)
  to_add_times <- matrix(NA, ncol = ncol(features), nrow = n_entries)

  colnames(to_add_ftrs) <- colnames(aligned$int_crosstab)
  colnames(to_add_times) <- colnames(aligned$rt_crosstab)

  sel <- seq_len(nrow(known_table))
  if (nrow(pairing) > 0) {
    sel <- sel[-(pairing[, 2])]
  }


  to_add_ftrs[, 1] <- to_add_times[, 1] <- known_table[sel, 6]
  to_add_ftrs[, 2] <- to_add_times[, 2] <- known_table[sel, 11]
  to_add_ftrs[, 3] <- to_add_times[, 3] <- known_table[sel, 9]
  to_add_ftrs[, 4] <- to_add_times[, 4] <- known_table[sel, 10]

  list(
    rt_crosstab = rbind(aligned$rt_crosstab, to_add_times),
    int_crosstab = rbind(aligned$int_crosstab, to_add_ftrs)
  )
}

augment_known_table <- function(
  aligned,
  known_table,
  match_tol_ppm,
  new_feature_min_count
) {
  pairing <- .merge_peaks(aligned, known_table, match_tol_ppm)
  rt_crosstab <- as.matrix(aligned$rt_crosstab)
  int_crosstab <- as.matrix(aligned$int_crosstab)

  for (i in seq_len(nrow(pairing))) {
    known_table[pairing[i, 2], ] <- peak.characterize(
      existing.row = known_table[pairing[i, 2], ],
      ftrs.row = int_crosstab[pairing[i, 1], ],
      rt.row = rt_crosstab[pairing[i, 1], ])
  }

  newly_found_ftrs <- which(!(seq_len(nrow(int_crosstab)) %in% pairing[, 1]))
  num_exp_found <- apply(int_crosstab != 0, 1, sum)

  for (i in newly_found_ftrs) {
    if (num_exp_found[i] >= new_feature_min_count) {
      row <- peak.characterize(
        existing.row = NA,
        ftrs.row = int_crosstab[i, ],
        rt.row = rt_crosstab[i, ])
      known_table <- rbind(known_table, row)
      pairing <- rbind(pairing, c(i, nrow(known_table)))
    }
  }

  list(pairing = pairing, known_table = known_table)
}

#' Runs features extraction in hybrid mode.
#' 
#' features extraction in hybrid mode.
#' 
#' @param filenames The CDF file names.
#' @param known_table Table of known chemicals.
#' @param min_exp A feature has to show up in at least this number of profiles to be included in the final result.
#' @param min_pres This is a parameter of the run filter, to be passed to the function proc.cdf().
#' @param min_run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param mz_tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. 
#'  This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy level. 
#'  Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline_correct After grouping the observations, the highest intensity in each group is found. If the highest is lower than 
#'  this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the 
#'  height of the noise groups.
#' @param baseline_correct_noise_percentile The perenctile of signal strength of those EIC that don't pass the run filter, to be used 
#'  as the baseline threshold of signal strength.
#' @param shape_model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian". When the 
#'  peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param BIC_factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1, models 
#'  with more peaks are penalized more.
#' @param peak_estim_method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param min_bandwidth The minimum bandwidth to use in the kernel smoother.
#' @param max_bandwidth The maximum bandwidth to use in the kernel smoother.
#' @param sd_cut A vector of two. Features with standard deviation outside the range defined by the two numbers are eliminated.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation and 
#'  the righ-standard deviation of the bi-Gaussian function used to fit the data.
#' @param component_eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a 
#'  proportion of intensities less than this value, the component will be ignored.
#' @param moment_power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param align_mz_tol The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the 
#'  tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z 
#'  value, becomes the cutoff level.
#' @param align_rt_tol The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for 
#'  the tolerance level based on the data.
#' @param max_align_mz_diff As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. 
#'  This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.
#' @param match_tol_ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param new_feature_min_count The number of profiles a new feature must be present for it to be added to the database.
#' @param recover_mz_range The m/z around the feature m/z to search for observations. The default value is NA, in which case 1.5 times 
#'  the m/z tolerance in the aligned object will be used.
#' @param recover_rt_range The retention time around the feature retention time to search for observations. The default value is NA, 
#'  in which case 0.5 times the retention time tolerance in the aligned object will be used.
#' @param use_observed_range If the value is TRUE, the actual range of the observed locations of the feature in all the spectra will be used.
#' @param recover_min_count Minimum number of raw data points to support a recovery.
#' @param intensity_weighted Whether to use intensity to weight mass density estimation.
#' @param cluster The number of CPU cores to be used
#' @export
#' @examples
#' hybrid(test_files, known_table, cluster = num_workers)
hybrid <- function(
  filenames,
  known_table,
  min_exp = 2,
  min_pres = 0.5,
  min_run = 12,
  mz_tol = 1e-05,
  baseline_correct = 0,
  baseline_correct_noise_percentile = 0.05,
  shape_model = "bi-Gaussian",
  BIC_factor = 2,
  peak_estim_method = "moment",
  min_bandwidth = NA,
  max_bandwidth = NA,
  sd_cut = c(0.01, 500),
  sigma_ratio_lim = c(0.01, 100),
  component_eliminate = 0.01,
  moment_power = 1,
  align_mz_tol = NA,
  align_rt_tol = NA,
  max_align_mz_diff = 0.01,
  match_tol_ppm = NA,
  new_feature_min_count = 2,
  recover_mz_range = NA,
  recover_rt_range = NA,
  use_observed_range = TRUE,
  recover_min_count = 3,
  intensity_weighted = FALSE,
  cluster = 4
) {
  if (!is(cluster, 'cluster')) {
    cluster <- parallel::makeCluster(cluster)
    on.exit(parallel::stopCluster(cluster))
  }

  # NOTE: side effect (doParallel has no functionality to clean up)
  doParallel::registerDoParallel(cluster)

  check_files(filenames)
  sample_names <- get_sample_name(filenames)

  message("**** feature extraction ****")
  extracted <- extract_features(
    cluster = cluster,
    filenames = filenames,
    min_pres = min_pres,
    min_run = min_run,
    mz_tol = mz_tol,
    baseline_correct = baseline_correct,
    baseline_correct_noise_percentile = baseline_correct_noise_percentile,
    intensity_weighted = intensity_weighted,
    min_bandwidth = min_bandwidth,
    max_bandwidth = max_bandwidth,
    sd_cut = sd_cut,
    sigma_ratio_lim = sigma_ratio_lim,
    shape_model = shape_model,
    peak_estim_method = peak_estim_method,
    component_eliminate = component_eliminate,
    moment_power = moment_power,
    BIC_factor = BIC_factor
  )

  message("**** time correction ****")
  corrected <- adjust.time(
    extracted_features = extracted,
    mz_tol_relative = align_mz_tol,
    rt_tol_relative = align_rt_tol,
    mz_max_diff = 10 * mz_tol,
    mz_tol_absolute = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** feature alignment ****")
  aligned <- align_features(
    sample_names = sample_names,
    features = corrected,
    min_occurrence = min_exp,
    mz_tol_relative = align_mz_tol,
    rt_tol_relative = align_rt_tol,
    mz_max_diff = 10 * mz_tol,
    mz_tol_absolute = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** augmenting with known peaks ****")
  merged <- augment_with_known_features(
    aligned = aligned,
    known_table = known_table,
    match_tol_ppm = match_tol_ppm
  )


  message("**** weaker signal recovery ****")
  recovered <- recover_weaker_signals(
    cluster = cluster,
    filenames = filenames,
    extracted_features = extracted,
    corrected_features = corrected,
    aligned_rt_crosstab = merged$rt_crosstab,
    aligned_int_crosstab = merged$int_crosstab,
    original_mz_tolerance = mz_tol,
    aligned_mz_tolerance = aligned$mz_tolerance,
    aligned_rt_tolerance = aligned$rt_tolerance,
    mz_range = recover_mz_range,
    rt_range = recover_rt_range,
    use_observed_range = use_observed_range,
    min_bandwidth = min_bandwidth,
    max_bandwidth = max_bandwidth,
    recover_min_count = recover_min_count
  )

  message("**** second round time correction ****")
  recovered_corrected <- adjust.time(
    extracted_features = recovered$extracted_features,
    mz_tol_relative = align_mz_tol,
    rt_tol_relative = align_rt_tol,
    mz_max_diff = 10 * mz_tol,
    mz_tol_absolute = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** second round feature alignment ****")
  recovered_aligned <- align_features(
    sample_names = sample_names,
    features = recovered_corrected,
    min_occurrence = min_exp,
    mz_tol_relative = align_mz_tol,
    rt_tol_relative = align_rt_tol,
    mz_max_diff = 10 * mz_tol,
    mz_tol_absolute = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** augmenting known table ****")
  augmented <- augment_known_table(
    aligned = recovered_aligned,
    known_table = known_table,
    match_tol_ppm = match_tol_ppm,
    new_feature_min_count = new_feature_min_count
  )

  aligned_feature_sample_table <- as_feature_sample_table(
    rt_crosstab = aligned$rt_crosstab,
    int_crosstab = aligned$int_crosstab
  )
  recovered_feature_sample_table <- as_feature_sample_table(
    rt_crosstab = recovered_aligned$rt_crosstab,
    int_crosstab = recovered_aligned$int_crosstab
  )

  list(
    extracted_features = recovered$extracted_features,
    corrected_features = recovered_corrected,
    aligned_feature_sample_table = aligned_feature_sample_table,
    recovered_feature_sample_table = recovered_feature_sample_table,
    aligned_mz_tolerance = as.numeric(recovered_aligned$mz_tolerance),
    aligned_rt_tolerance = as.numeric(recovered_aligned$rt_tolerance),
    updated_known_table = as.data.frame(augmented$known_table),
    features_known_table_pairing = as.data.frame(augmented$pairing)
  )
}
