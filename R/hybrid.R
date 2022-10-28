#' @import parallel doParallel
NULL
#> NULL

#' Match peaks from sample table to already known peaks via similar m/z and rt.
#' @param aligned A list object with three tibble tables: metadata, intensity, and rt.
#' @param known_table A table of known/previously detected peaks.
#' @param match_tol_ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the 
#'  tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z 
#'  value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for 
#'  the tolerance level based on the data.
#' @return n x 2 matrix containing sample features-known features pairs.
match_peaks <- function(aligned,
  known_table,
  match_tol_ppm,
  mz_tol_relative,
  rt_tol_relative) {

  if (is.na(match_tol_ppm)) {
    match_tol_ppm <- mz_tol_relative * 1e+06
  }

  mass_matched_pos <- find_mz_match(aligned$metadata[['mz']],
    known_table['m.z'],
    match_tol_ppm)

  known_assigned <- rep(0, nrow(known_table))
  new_assigned <- rep(0, nrow(aligned$metadata))
  pairing <- matrix(0, nrow = 0, ncol = 2)

  for (i in mass_matched_pos) {
    if (new_assigned[i] != 0) {
      next
    }
    # find all potentially related known/newly found peaks
    prev_sel_new <- i
    threshold <- aligned$metadata[[i, 'mz']] * mz_tol_relative

    sel_known <- which(abs(known_table[['m.z']] - aligned$metadata[[i, 'mz']]) < threshold)
    sel_new <- c()
    for (m in seq_along(sel_known)) {
      distance <- abs(aligned$metadata[['mz']] - known_table[[sel_known[m], 'm.z']])
      sel_new <- unique(c(sel_new, which(distance < threshold)))
    }

    while (length(sel_new) > length(prev_sel_new)) {
      prev_sel_new <- sel_new

      sel_known <- NULL
      for (m in seq_along(sel_new)) {
        distance <- abs(known_table[['m.z']] - aligned$metadata[[sel_new[m], 'mz']])
        sel_known <- unique(c(sel_known, which(distance < threshold)))
      }

      sel_new <- NULL
      for (m in seq_along(sel_known)) {
        distance <- abs(aligned$metadata[['mz']] - known_table[[sel_known[m], 'm.z']])
        sel_new <- unique(c(sel_new, which(distance < threshold)))
      }
    }

    time_matched <- mass_matched <- matrix(
      data = 0,
      nrow = length(sel_known),
      ncol = length(sel_new))

    for (k in seq_along(sel_known)) {
      time_matched[k,] <- abs(aligned$metadata$rt[sel_new] - known_table[[sel_known[k], 'RT_mean']])
      mass_matched[k,] <- abs(aligned$metadata$mz[sel_new] - known_table[[sel_known[k], 'm.z']])
    }

    mass_matched <- mass_matched/median(known_table[sel_known, 'm.z'])
    time_matched[mass_matched <= match_tol_ppm * 1e-06] <- 1e+10

    time_matched[is.na(time_matched)] <- rt_tol_relative / 2
    both_matched <- find.match(time_matched, rt_tol_relative / 2)

    for (m in seq_along(sel_new)) {
      k <- which(both_matched[, m] == 1)

      if (length(k) == 1 && known_assigned[sel_known[k]] == 0) {
        new_assigned[sel_new[m]] <- 1
        known_assigned[sel_known[k]] <- 1
        pairing <- rbind(pairing, c(sel_new[m], sel_known[k]))
      }
    }
  }
  colnames(pairing) <- c('new', 'known')
  return(pairing)
}


#' A wrapper function to join knowledge from aligned features and known table.
#' 
#' @param features A list object with three tibble tables: metadata, intensity, and rt.
#' @param known_table A table of known/previously detected peaks.
#' @param match_tol_ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the 
#'  tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z 
#'  value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for 
#'  the tolerance level based on the data.
#' @param from_features_to_known Determines direction of joining; if TRUE, aligned features are joined to known table, vice verse if it is FALSE.
#' @param new_feature_min_count The number of profiles a new feature must be present for it to be added to the database.
#' @return Enriched aligned table or known features.
#' @import dplyr
#' @export
merge_features_and_known_table <- function(
  features,
  known_table,
  match_tol_ppm,
  mz_tol_relative,
  rt_tol_relative,
  from_features_to_known_table = TRUE,
  new_feature_min_count = NA) {
    if (from_features_to_known_table) {
        return(augment_known_table(features,
                                   known_table,
                                   match_tol_ppm,
                                   mz_tol_relative,
                                   rt_tol_relative,
                                   new_feature_min_count)
               )
    } else {
        return(enrich_table_by_known_features(features,
                                              known_table,
                                              match_tol_ppm,
                                              mz_tol_relative,
                                              rt_tol_relative)
               )
    }
}


#' Add entries from the known features table to the aligned table.
#' 
#' @param aligned A list object with three tibble tables: metadata, intensity, and rt.
#' @param known_table A table of known/previously detected peaks.
#' @param match_tol_ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the 
#'  tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z 
#'  value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for 
#'  the tolerance level based on the data.
#' @return Aligned table with known features.
#' @import dplyr
enrich_table_by_known_features <- function(
  aligned,
  known_table,
  match_tol_ppm,
  mz_tol_relative,
  rt_tol_relative
  ) {
  pairing <- match_peaks(aligned, known_table, match_tol_ppm, mz_tol_relative, rt_tol_relative)

  known_table <- tibble(known_table)[-pairing[,'known'], ]
  metadata <- select(known_table, c('m.z', 'mz_min', 'mz_max', 'RT_mean', 'RT_min', 'RT_max'))
  colnames(metadata) <- c('mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax')

  new_features_num <- nrow(metadata)
  samples_num <- ncol(select(aligned$intensity, -id))

  aligned$metadata <- bind_rows(select(aligned$metadata, -id), metadata) |>
    rowid_to_column('id')

  rt <- data.frame(matrix(data = NA, nrow = new_features_num, ncol = samples_num))
  colnames(rt) <- colnames(select(aligned$rt, -id))
  aligned$rt <- bind_rows(select(aligned$rt, -id), rt) |>
    rowid_to_column('id')

  intensity <- data.frame(matrix(data = 0, nrow = new_features_num, ncol = samples_num))
  colnames(intensity) <- colnames(select(aligned$intensity, -id))
  aligned$intensity <- bind_rows(select(aligned$intensity, -id), intensity) |>
    rowid_to_column('id')

  return(aligned)
}

#' Add newly detected aligned features to a known features table.
#' 
#' @param aligned A list object with three tibble tables: metadata, intensity, and rt.
#' @param known_table A table of known/previously detected peaks.
#' @param match_tol_ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the 
#'  tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z 
#'  value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for 
#'  the tolerance level based on the data.
#' @param new_feature_min_count The number of profiles a new feature must be present for it to be added to the database.
#' @return Known table with novel features.
augment_known_table <- function(
  aligned,
  known_table,
  match_tol_ppm,
  mz_tol_relative,
  rt_tol_relative,
  new_feature_min_count
) {
  pairing <- match_peaks(aligned, known_table, match_tol_ppm, mz_tol_relative, rt_tol_relative)

  for (i in seq_len(nrow(pairing))) {
    known_table[pairing[i, 'known'], ] <- peak_characterize(
      existing_row = known_table[pairing[i, 'known'], ],
      metadata_row = aligned$metadata[pairing[i, 'new'], ],
      ftrs_row = aligned$intensity[pairing[i, 'new'], ],
      rt_row = aligned$rt[pairing[i, 'new'], ])
  }

  newly_found_ftrs <- which(!(seq_len(nrow(aligned$metadata)) %in% pairing[, 'new']))
  num_exp_found <- apply(aligned$intensity != 0, 1, sum)

  for (i in newly_found_ftrs) {
    if (num_exp_found[i] >= new_feature_min_count) {
      row <- peak_characterize(
        existing_row = NA,
        metadata_row = aligned$metadata[i, ],
        ftrs_row = aligned$intensity[i, ],
        rt_row = aligned$rt[i, ])
      known_table <- dplyr::bind_rows(known_table, row)
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
  register_functions_to_cluster(cluster)

  check_files(filenames)
  sample_names <- get_sample_name(filenames)
  number_of_samples <- length(sample_names)
  
  message("**** feature extraction ****")
  profiles <- snow::parLapply(cluster, filenames, function(filename) {
      proc.cdf(
          filename = filename,
          min_presence = min_pres,
          min_elution_length = min_run,
          mz_tol = mz_tol,
          baseline_correct = baseline_correct,
          baseline_correct_noise_percentile = baseline_correct_noise_percentile,
          intensity_weighted = intensity_weighted,
          do.plot = FALSE,
          cache = FALSE
      )
  })
  
  extracted <- snow::parLapply(cluster, profiles, function(profile) {
      prof.to.features(
          profile = profile,
          min.bw = min_bandwidth,
          max.bw = max_bandwidth,
          sd.cut = sd_cut,
          sigma.ratio.lim = sigma_ratio_lim,
          shape.model = shape_model,
          estim.method = peak_estim_method,
          component.eliminate = component_eliminate,
          power = moment_power,
          BIC.factor = BIC_factor,
          do.plot = FALSE
      )
  })

 message("**** computing clusters ****")
  extracted_clusters <- compute_clusters(
    feature_tables = extracted,
    mz_tol_relative = align_mz_tol,
    mz_tol_absolute = max_align_mz_diff,
    mz_max_diff = 10 * mz_tol,
    rt_tol_relative = align_rt_tol,
    sample_names = sample_names
  )

  message("**** computing template ****")
  template_features <- compute_template(extracted_clusters$feature_tables)


  message("**** time correction ****")
  corrected <- foreach::foreach(this.feature = extracted_clusters$feature_tables) %do% correct_time(
    this.feature,
    template_features,
    extracted_clusters$mz_tol_relative,
    extracted_clusters$rt_tol_relative
  )

  message("**** computing clusters ****")
  adjusted_clusters <- compute_clusters(
    feature_tables = corrected,
    mz_tol_relative = extracted_clusters$mz_tol_relative,
    mz_tol_absolute = extracted_clusters$rt_tol_relative,
    mz_max_diff = 10 * mz_tol,
    rt_tol_relative = align_rt_tol
  )

  message("**** feature alignment ****")
  aligned <- create_aligned_feature_table(
      dplyr::bind_rows(adjusted_clusters$feature_tables),
      min_exp,
      sample_names,
      adjusted_clusters$rt_tol_relative,
      adjusted_clusters$mz_tol_relative
  )

  message("**** augmenting with known peaks ****")
  merged <- merge_features_and_known_table(
      features = aligned,
      known_table = known_table,
      match_tol_ppm = match_tol_ppm,
      mz_tol_relative = adjusted_clusters$mz_tol_relative,
      rt_tol_relative = adjusted_clusters$rt_tol_relative,
      from_features_to_known_table = FALSE
  )

  message("**** weaker signal recovery ****")
  recovered <- lapply(seq_along(filenames), function(i) {
    recover.weaker(
      filename = filenames[[i]],
      sample_name = as.character(i),
      extracted_features = extracted[[i]],
      adjusted_features = corrected[[i]],
      metadata_table = merged$metadata,
      rt_table = merged$rt,
      intensity_table = merged$intensity,
      orig.tol = mz_tol,
      align.mz.tol = extracted_clusters$mz_tol_relative,
      align.rt.tol = extracted_clusters$rt_tol_relative,
      recover_mz_range = recover_mz_range,
      recover_rt_range = recover_rt_range,
      use.observed.range = use_observed_range,
      bandwidth = 0.5,
      min.bw = min_bandwidth,
      max.bw = max_bandwidth,
      recover.min.count = recover_min_count,
      intensity.weighted = intensity_weighted
    )
  })

  recovered_adjusted <- lapply(recovered, function(x) x$adjusted_features)

  message("**** third time computing clusters ****")
  recovered_clusters <- compute_clusters(
    feature_tables = recovered_adjusted,
    mz_tol_relative = align_mz_tol,
    mz_tol_absolute = max_align_mz_diff,
    mz_max_diff = 10 * mz_tol,
    rt_tol_relative = align_rt_tol
  )

  message("**** computing template ****")
  template_features <- compute_template(recovered_clusters$feature_tables)


  message("**** second time correction ****")
  corrected <- foreach::foreach(this.feature = recovered_clusters$feature_tables) %do% correct_time(
    this.feature,
    template_features,
    recovered_clusters$mz_tol_relative,
    recovered_clusters$rt_tol_relative
  )

  message("**** fourth computing clusters ****")
  adjusted_clusters <- compute_clusters(
    feature_tables = corrected,
    mz_tol_relative = recovered_clusters$mz_tol_relative,
    mz_tol_absolute = recovered_clusters$rt_tol_relative,
    mz_max_diff = 10 * mz_tol,
    rt_tol_relative = align_rt_tol
  )

  message("**** second feature alignment ****")
  recovered_aligned <- create_aligned_feature_table(
      dplyr::bind_rows(adjusted_clusters$feature_tables),
      min_exp,
      sample_names,
      adjusted_clusters$rt_tol_relative,
      adjusted_clusters$mz_tol_relative
  )

  message("**** augmenting known table ****")
  augmented <- merge_features_and_known_table(
    features = recovered_aligned,
    known_table = known_table,
    match_tol_ppm = match_tol_ppm,
    mz_tol_relative = adjusted_clusters$mz_tol_relative,
    rt_tol_relative = adjusted_clusters$rt_tol_relative,
    from_features_to_known_table = TRUE,
    new_feature_min_count = new_feature_min_count
  )

  aligned_feature_sample_table <- as_feature_sample_table(
    metadata = aligned$metadata,
    rt_crosstab = aligned$rt,
    int_crosstab = aligned$intensity
  )
  recovered_feature_sample_table <- as_feature_sample_table(
    metadata = recovered_aligned$metadata,
    rt_crosstab = recovered_aligned$rt,
    int_crosstab = recovered_aligned$intensity
  )

  list(
    extracted_features = recovered$extracted_features,
    corrected_features = corrected,
    aligned_feature_sample_table = aligned_feature_sample_table,
    recovered_feature_sample_table = recovered_feature_sample_table,
    aligned_mz_tolerance = as.numeric(recovered_aligned$mz_tolerance),
    aligned_rt_tolerance = as.numeric(recovered_aligned$rt_tolerance),
    updated_known_table = as.data.frame(augmented$known_table),
    features_known_table_pairing = as.data.frame(augmented$pairing)
  )
}
