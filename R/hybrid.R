.merge_peaks <- function(aligned, known_table, match_tol_ppm) {
  if (is.na(match_tol_ppm)) {
    match_tol_ppm <- aligned$mz_tolerance * 1e+06
  }

  features <- as.matrix(aligned$int_crosstab)
  known_mz <- known_table[, 6]
  known_rt <- known_table[, 11]

  mass_d2 <- mass.match(features[, 1], known_mz, match_tol_ppm)
  mass_matched_pos <- which(mass_d2 > 0)

  known_assigned <- rep(0, nrow(known_table))
  new_assigned <- rep(0, nrow(features))
  pairing <- matrix(0, ncol = 2, nrow = 0)
  colnames(pairing) <- c("new", "known")

  for (i in mass_matched_pos) {
    if (new_assigned[i] == 0) {
      # find all potentially related known/newly found peaks
      prev_sel_new <- i
      threshold <- features[i, 1] * match_tol_ppm / 1e+06

      sel_known <- which(abs(known_mz - features[i, 1]) < threshold)
      sel_new <- NULL
      for (m in seq_along(sel_known)) {
        distance <- abs(features[, 1] - known_mz[sel_known[m]])
        sel_new <- c(sel_new, which(distance < threshold))
      }
      sel_known <- unique(sel_known)
      sel_new <- unique(sel_new)

      while (length(sel_new) > length(prev_sel_new)) {
        prev_sel_new <- sel_new

        sel_known <- NULL
        for (m in seq_along(sel_new)) {
          distance <- abs(known_mz - features[sel_new[m], 1])
          sel_known <- c(sel_known, which(distance < threshold))
        }

        sel_new <- NULL
        for (m in seq_along(sel_known)) {
          distance <- abs(features[, 1] - known_mz[sel_known[m]])
          sel_new <- c(sel_new, which(distance < threshold))
        }

        sel_known <- unique(sel_known)
        sel_new <- unique(sel_new)
      }

      time_matched <- mass_matched <-
        matrix(data = 0, nrow = length(sel_known), ncol = length(sel_new))

      for (k in seq_along(sel_known)) {
        time_matched[k, ] <- abs(features[sel_new, 2] - known_rt[sel_known[k]])
        mass_matched[k, ] <- abs(features[sel_new, 1] - known_mz[sel_known[k]])
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
      chr.row = rt_crosstab[pairing[i, 1], ])
  }

  newly_found_ftrs <- which(!(seq_len(nrow(int_crosstab)) %in% pairing[, 1]))
  num_exp_found <- apply(int_crosstab != 0, 1, sum)

  for (i in newly_found_ftrs) {
    if (num_exp_found[i] >= new_feature_min_count) {
      row <- peak.characterize(
        existing.row = NA,
        ftrs.row = int_crosstab[i, ],
        chr.row = rt_crosstab[i, ])
      known_table <- rbind(known_table, row)
      pairing <- rbind(pairing, c(i, nrow(known_table)))
    }
  }

  list(pairing = pairing, known_table = known_table)
}

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
  align_chr_tol = NA,
  max_align_mz_diff = 0.01,
  match_tol_ppm = NA,
  new_feature_min_count = 2,
  recover_mz_range = NA,
  recover_chr_range = NA,
  use_observed_range = TRUE,
  recover_min_count = 3,
  intensity_weighted = FALSE,
  cluster = detectCores()
) {
  if (!is(cluster, 'cluster')) {
    cluster <- makeCluster(cluster)
    on.exit(parallel::stopCluster(cluster))
  }

  # NOTE: side effect (doParallel has no functionality to clean up)
  doParallel::registerDoParallel(cluster)

  check_files(filenames)
  # further processing requires sorted file list according to the acquisition
  # order (assumes that the filenames contain acquisition number)
  filenames <- sort(unlist(filenames))
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
    features = extracted,
    mz.tol = align_mz_tol,
    chr.tol = align_chr_tol,
    find.tol.max.d = 10 * mz_tol,
    max.align.mz.diff = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** feature alignemnt ****")
  aligned <- align_features(
    sample_names = sample_names,
    features = corrected,
    min.exp = min_exp,
    mz.tol = align_mz_tol,
    chr.tol = align_chr_tol,
    find.tol.max.d = 10 * mz_tol,
    max.align.mz.diff = max_align_mz_diff,
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
    rt_range = recover_chr_range,
    use_observed_range = use_observed_range,
    min_bandwidth = min_bandwidth,
    max_bandwidth = max_bandwidth,
    recover_min_count = recover_min_count
  )

  message("**** second round time correction ****")
  recovered_corrected <- adjust.time(
    features = recovered$extracted_features,
    mz.tol = align_mz_tol,
    chr.tol = align_chr_tol,
    find.tol.max.d = 10 * mz_tol,
    max.align.mz.diff = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** second round feature alignment ****")
  recovered_aligned <- align_features(
    sample_names = sample_names,
    features = recovered_corrected,
    min.exp = min_exp,
    mz.tol = align_mz_tol,
    chr.tol = align_chr_tol,
    find.tol.max.d = 10 * mz_tol,
    max.align.mz.diff = max_align_mz_diff,
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
    aligned_mz_toletance = as.numeric(recovered_aligned$mz_tolerance),
    aligned_rt_tolerance = as.numeric(recovered_aligned$rt_tolerance),
    updated_known_table = as.data.frame(augmented$known_table),
    features_known_table_pairing = as.data.frame(augmented$pairing)
  )
}
