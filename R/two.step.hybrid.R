filter_features_by_presence <- function(
  feature_table, 
  metadata, 
  batches_idx,
  within_batch_threshold,
  across_batch_threshold) {
  across_batch_presence <- data.frame(
    matrix(
      nrow = nrow(feature_table), 
      ncol = length(batches_idx)))
  for (batch_id in batches_idx) {
    samples <- filter(metadata, batch == batch_id)$sample_name
    intensities <- dplyr::select(feature_table, samples)
    presence <- intensities > 0
    above_threshold <- rowSums(presence) / length(samples) >= within_batch_threshold
    across_batch_presence[,batch_id] <- above_threshold
  }
  above_threshold <- rowSums(across_batch_presence) / length(batches_idx) >= across_batch_threshold
  return(feature_table[above_threshold, ])
}

as_wide_aligned_table <- function(aligned) {
  mz_scale_table <- aligned$rt_crosstab[, c("mz", "rt", "mz_min", "mz_max")]
  aligned <- as_feature_sample_table(
    rt_crosstab = aligned$rt_crosstab,
    int_crosstab = aligned$int_crosstab
  )
  aligned <- long_to_wide_feature_table(aligned)
  aligned <- dplyr::inner_join(aligned, mz_scale_table, by = c("mz", "rt"))
  return(aligned)
}

pivot_feature_values <- function(feature_table, variable) {
  extended_variable <- paste0("sample_", variable)
  values <- dplyr::select(feature_table, mz, rt, sample, !!sym(extended_variable))
  values <- tidyr::pivot_wider(values, names_from = sample, values_from = !!sym(extended_variable))
  variable_colnames <- colnames(values)[3:ncol(values)]
  variable_colnames <- paste0(variable_colnames, "_", variable)
  colnames(values)[3:ncol(values)] <- variable_colnames
  return(values)
}

long_to_wide_feature_table <- function(feature_table) {
  sample_rts <- pivot_feature_values(feature_table, "rt")
  sample_intensities <- pivot_feature_values(feature_table, "intensity")
  feature_table <- dplyr::select(feature_table, mz, rt) %>%
    dplyr::distinct(mz, rt) %>%
    dplyr::inner_join(sample_rts, by = c("mz", "rt")) %>%
    dplyr::inner_join(sample_intensities, by = c("mz", "rt"))
}

readjust_times <- function(within_batch, between_batch) {
  within_batch_recovered <- long_to_wide_feature_table(
    within_batch$recovered_feature_sample_table)
  between_batch_rts <- between_batch$rt
  for (j in 1:length(within_batch$corrected_features)) {
    for (i in 1:nrow(within_batch$corrected_features[[j]])) {
      diff.time <- abs(
        within_batch_recovered$rt -
        within_batch$corrected_features[[j]][i, "pos"])
      min_idx <- which(diff.time == min(diff.time))[1]
      within_batch$corrected_features[[j]][i, "pos"] <- between_batch_rts[min_idx]
    }
  }
  return(within_batch$corrected_features)
}

compute_intensity_medians <- function(feature_table) {
  stopifnot("sample_intensity" %in% colnames(feature_table))
  feature_table <- dplyr::group_by(feature_table, feature, mz, rt) %>%
    dplyr::mutate(intensity = median(sample_intensity)) %>%
    dplyr::ungroup()
  return(feature_table)
}

bind_batch_label_column <- function(filenames, metadata) {
  stopifnot(nrow(metadata) == length(filenames))

  filenames <- as_tibble(filenames)
  colnames(filenames)[1] <- "filename"
  filenames <- mutate(filenames, sample_name = get_sample_name(filename))
  filenames <- inner_join(filenames, metadata, by = "sample_name")
  return (dplyr::select(filenames, -sample_name))
}

two.step.hybrid <- function(
  filenames,
  metadata,
  min.within.batch.prop.detect = 0.1,
  min.within.batch.prop.report = 0.5,
  min.batch.prop = 0.5,
  batch.align.mz.tol = 1e-5,
  batch.align.chr.tol = 50,
  known.table = NA,
  cluster = 4,
  min.pres = 0.5,
  min.run = 12,
  mz.tol = 1e-5,
  baseline.correct.noise.percentile = 0.05,
  shape.model = "bi-Gaussian",
  baseline.correct = 0,
  peak.estim.method = "moment",
  min.bw = NA,
  max.bw = NA,
  sd.cut = c(0.1, 100),
  sigma.ratio.lim = c(0.05, 20),
  component.eliminate = 0.01,
  moment.power = 2,
  align.mz.tol = NA,
  align.chr.tol = NA,
  max.align.mz.diff = 0.01,
  pre.process = FALSE,
  recover.mz.range = NA,
  recover.chr.range = NA,
  use.observed.range = TRUE,
  match.tol.ppm = NA,
  new.feature.min.count = 2,
  recover.min.count = 3,
  intensity.weighted = FALSE,
  BIC.factor = 2) 
  {

  metadata <- read.table(metadata, sep=",", header=TRUE)
  filenames_batchwise <- bind_batch_label_column(filenames, metadata)
  batches_idx <- unique(metadata$batch)

  batchwise <- new("list")
  message("**** processing ", length(batches_idx), " batches separately ****")
  for (batch_id in batches_idx) {
    filenames_batch <- dplyr::filter(filenames_batchwise, batch == batch_id)$filename
    batchwise_features <- hybrid(
      filenames = filenames_batch,
      known_table = known.table,
      min_exp = ceiling(min.within.batch.prop.detect * length(filenames_batchwise)),
      min_pres = min.pres,
      min_run = min.run,
      mz_tol = mz.tol,
      baseline_correct = baseline.correct,
      baseline_correct_noise_percentile = baseline.correct.noise.percentile,
      shape_model = shape.model,
      BIC_factor = BIC.factor,
      peak_estim_method = peak.estim.method,
      min_bandwidth = min.bw,
      max_bandwidth = max.bw,
      sd_cut = sd.cut,
      sigma_ratio_lim = sigma.ratio.lim,
      component_eliminate = component.eliminate,
      moment_power = moment.power,
      align_mz_tol = align.mz.tol,
      align_chr_tol = align.chr.tol,
      max_align_mz_diff = max.align.mz.diff,
      match_tol_ppm = match.tol.ppm,
      new_feature_min_count = new.feature.min.count,
      recover_mz_range = recover.mz.range,
      recover_chr_range = recover.chr.range,
      use_observed_range = use.observed.range,
      recover_min_count = recover.min.count,
      intensity_weighted = intensity.weighted,
      cluster = cluster
    )

    batchwise[[batch_id]] <- batchwise_features
  }
  step_one_features <- list()
  for (batch_id in batches_idx) {
    step_one_features[[batch_id]] <- compute_intensity_medians(
      batchwise[[batch_id]]$recovered_feature_sample_table) %>%
      dplyr::select(-feature)
  }

  cl <- makeCluster(cluster)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(recetox.aplcms))

  message("*** aligning time ***")
  corrected <- adjust.time(step_one_features,
    mz.tol = batch.align.mz.tol,
    chr.tol = batch.align.chr.tol,
    find.tol.max.d = 10 * mz.tol,
    max.align.mz.diff = max.align.mz.diff,
    rt_colname = "rt")

  message("*** aligning features ***")
  aligned <- align_features(
    sample_names = paste0("batch_", batches_idx),
    features = corrected,
    min.exp = ceiling(min.batch.prop * length(batches_idx)),
    mz.tol = batch.align.mz.tol,
    chr.tol = batch.align.chr.tol,
    find.tol.max.d = 10 * mz.tol,
    max.align.mz.diff = max.align.mz.diff,
    rt_colname = "rt")

  aligned <- as_wide_aligned_table(aligned)

  stopCluster(cl)

  message("*** recovering features across batches ***")

  cl <- makeCluster(cluster)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(recetox.aplcms))


  recovered <- tibble(mz = numeric(),
    rt = numeric(),
    mz_min = numeric(),
    mz_max = numeric())

  for (batch_id in batches_idx)
  {
    this.fake <- step_one_features[[batch_id]]
    this.fake.medians <- distinct(this.fake, mz, rt, intensity)$intensity
    this.fake <- long_to_wide_feature_table(this.fake)


    # adjusting the time (already within batch adjusted)
    this.features <- readjust_times(batchwise[[batch_id]], corrected[[batch_id]])

    aligned_intensities <- dplyr::select(aligned, contains("_intensity"))
    batchwise_intensities <- dplyr::select(this.fake, contains("_intensity"))

    this.fake.time <- dplyr::select(this.fake, contains("_rt"))
    this.pk.time <- this.aligned <- matrix(0, nrow = nrow(aligned), ncol = ncol(batchwise_intensities))

    for (sample in 1:nrow(aligned)) {
      if (aligned_intensities[sample, batch_id] != 0) {
        idx <- which(between(this.fake$mz, aligned[sample, "mz_min"], aligned[sample, "mz_max"]) 
          & abs(this.fake.medians - aligned_intensities[sample, batch_id]) < 1)
        if (length(idx) < 1) {
          idx <- which(between(this.fake$mz, aligned[sample, "mz_min"], aligned[sample, "mz_max"]))
        }
        if (length(idx) < 1) {
          message("Warning: batch ", batch_id, " sample ", sample, " has matching issue.")
        } else {
          this.aligned[sample, ] <- apply(batchwise_intensities[idx, ], 2, sum)
          this.pk.time[sample, ] <- apply(this.fake.time[idx, ], 2, median)
        }
      } else {
        ### go into individual feature tables to find a match
        recaptured <- rep(0, ncol(this.aligned))
        recaptured.time <- rep(NA, ncol(this.aligned))

        for (j in 1:length(this.features)) {
          diff.mz <- abs(this.features[[j]][, "mz"] - aligned[sample, "mz"])
          diff.time <- abs(this.features[[j]][, "pos"] - aligned[sample, "rt"])
          idx <- which(diff.mz < aligned[sample, "mz"] * batch.align.mz.tol & diff.time <= batch.align.chr.tol)
          
          if (length(idx) > 0) {
            idx <- idx[which(diff.time[idx] == min(diff.time[idx]))[1]]
            recaptured[j] <- this.features[[j]][idx, "area"]
            recaptured.time[j] <- this.features[[j]][idx, "rt"]
          }
        }
        this.aligned[sample, ] <- recaptured
        this.pk.time[sample, ] <- recaptured.time
      }
    }

    colnames(this.aligned) <- extract_pattern_colnames(this.fake, "_intensity")
    colnames(this.aligned) <- stringr::str_remove_all(colnames(this.aligned), "_intensity")
    colnames(this.pk.time) <- extract_pattern_colnames(this.fake, "_rt")
    colnames(this.pk.time) <- stringr::str_remove_all(colnames(this.pk.time), "_rt")


    aligned_features <- dplyr::select(aligned, mz, rt, mz_min, mz_max)   
    aligned_int_crosstab <- dplyr::bind_cols(aligned_features, this.aligned)
    aligned_rt_crosstab <- dplyr::bind_cols(aligned_features, this.pk.time)
    recovered_batchwise <- recover_weaker_signals(
      cluster = cl,
      filenames = filter(filenames_batchwise, batch == batch_id)$filename,
      extracted_features = batchwise[[batch_id]]$extracted_features,
      corrected_features = batchwise[[batch_id]]$corrected_features,
      aligned_rt_crosstab = aligned_rt_crosstab,
      aligned_int_crosstab = aligned_int_crosstab,
      original_mz_tolerance = mz.tol,
      aligned_mz_tolerance = batch.align.mz.tol,
      aligned_rt_tolerance = batch.align.chr.tol,
      mz_range = recover.mz.range,
      rt_range = recover.chr.range,
      use_observed_range = use.observed.range,
      min_bandwidth = min.bw,
      max_bandwidth = max.bw,
      recover_min_count = recover.min.count
    )

    recovered_batchwise <- as_wide_aligned_table(recovered_batchwise)

    recovered <- dplyr::full_join(recovered, recovered_batchwise, by = c("mz", "rt", "mz_min", "mz_max"))
  }

  recovered <- dplyr::select(recovered, mz, rt, mz_min, mz_max, contains("_intensity"))
  stopCluster(cl)

  colnames(recovered) <- stringr::str_remove_all(colnames(recovered), "_intensity")
  aligned_filtered <- filter_features_by_presence(
    feature_table = recovered, 
    metadata = metadata, 
    batches_idx = batches_idx, 
    within_batch_threshold = min.within.batch.prop.report, 
    across_batch_threshold = min.batch.prop)

  features <- new("list")
  features$batchwise_results <- batchwise
  features$all_features <- aligned
  features$final_features <- recovered
  return(features)
}
