filter_features_by_presence <- function(feature_table,
                                        metadata,
                                        batches_idx,
                                        within_batch_threshold,
                                        across_batch_threshold) {
  across_batch_presence <- data.frame(
    matrix(
      nrow = nrow(feature_table),
      ncol = length(batches_idx)
    )
  )
  for (batch_id in batches_idx) {
    samples <- filter(metadata, batch == batch_id)$sample_name
    intensities <- dplyr::select(feature_table, contains("_intensity"))
    presence <- intensities > 0
    above_threshold <- rowSums(presence) / length(samples) >= within_batch_threshold
    across_batch_presence[, batch_id] <- above_threshold
  }
  above_threshold <- rowSums(across_batch_presence) / length(batches_idx) >= across_batch_threshold
  return(feature_table[above_threshold, ])
}

readjust_times <- function(within_batch, between_batch) {
  within_batch_recovered <- long_to_wide_feature_table(
    within_batch$recovered_feature_sample_table
  )
  between_batch_rts <- between_batch$rt
  for (j in 1:length(within_batch$corrected_features)) {
    for (i in 1:nrow(within_batch$corrected_features[[j]])) {
      diff.time <- abs(
        within_batch_recovered$rt -
          within_batch$corrected_features[[j]][i, "pos"]
      )
      min_idx <- which(diff.time == min(diff.time))[1]
      within_batch$corrected_features[[j]][i, "pos"] <- between_batch_rts[min_idx]
    }
  }
  return(within_batch$corrected_features)
}

compute_intensity_medians <- function(feature_table) {
  stopifnot("sample_intensity" %in% colnames(feature_table))
  feature_table <- dplyr::group_by(feature_table, mz, rt) %>%
    dplyr::mutate(median_intensity = median(sample_intensity)) %>%
    dplyr::ungroup()
  return(feature_table)
}

bind_batch_label_column <- function(filenames, metadata) {
  stopifnot(nrow(metadata) == length(filenames))

  filenames <- as_tibble(filenames)
  colnames(filenames)[1] <- "filename"
  filenames <- mutate(filenames, sample_name = get_sample_name(filename))
  filenames <- inner_join(filenames, metadata, by = "sample_name")
  return(dplyr::select(filenames, -sample_name))
}

feature_recovery <- function(cluster,
                             step_one_features,
                             batchwise,
                             filenames_batchwise,
                             corrected,
                             aligned,
                             batches_idx,
                             mz.tol,
                             batch.align.mz.tol,
                             batch.align.chr.tol,
                             recover.mz.range,
                             recover.chr.range,
                             use.observed.range,
                             min.bw,
                             max.bw,
                             recover.min.count) {
  recovered <- tibble(
    mz = numeric(),
    rt = numeric(),
    mz_min = numeric(),
    mz_max = numeric()
  )

  for (batch_id in batches_idx)
  {
    this.fake <- long_to_wide_feature_table(step_one_features[[batch_id]])
    this.fake.medians <- apply(dplyr::select(this.fake, contains("_intensity")), 1, median)

    # adjusting the time (already within batch adjusted)
    this.features <- readjust_times(batchwise[[batch_id]], corrected[[batch_id]])
    aligned_intensities <- dplyr::select(aligned, contains("_intensity"))
    batchwise_intensities <- dplyr::select(this.fake, contains("_intensity"))

    this.fake.time <- dplyr::select(this.fake, contains("_rt"))
    this.pk.time <- this.aligned <- matrix(0, nrow = nrow(aligned), ncol = ncol(batchwise_intensities))

    for (sample in 1:nrow(aligned)) {
      if (aligned_intensities[sample, batch_id] != 0) {
        idx <- which(between(this.fake$mz, aligned[sample, "mz_min"], aligned[sample, "mz_max"]) &
          abs(this.fake.medians - aligned_intensities[sample, batch_id]) < 1)
        if (length(idx) < 1) {
          idx <- which(between(this.fake$mz, aligned[sample, "mz_min"], aligned[sample, "mz_max"]))
        }
        if (length(idx) < 1) {
          message("warning: batch ", batch_id, " sample ", sample, " has matching issue")
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
            recaptured.time[j] <- this.features[[j]][idx, "pos"]
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
    aligned_int_crosstab <- dplyr::bind_cols(aligned_features, as_tibble(this.aligned))
    aligned_rt_crosstab <- dplyr::bind_cols(aligned_features, as_tibble(this.pk.time))
    recovered_batchwise <- recover_weaker_signals(
      cluster = cluster,
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
  recovered <- dplyr::select(recovered, mz, rt, mz_min, mz_max, contains("_rt"), contains("_intensity"))

  return(recovered)
}

semisup_to_hybrid_adapter <- function(batchwise, batches_idx) {
  for (batch in batches_idx) {
    final.ftrs <- as_tibble(batchwise[[batch]]$final.ftrs)
    final.times <- as_tibble(batchwise[[batch]]$final.times)

    mz_pattern <- c("mz.min", "mz.max")
    mz_replacement <- c("mz_min", "mz_max")

    colnames(final.ftrs) <- stringr::str_replace_all(
      colnames(final.ftrs),
      mz_pattern,
      mz_replacement)

    colnames(final.times) <- stringr::str_replace_all(
      colnames(final.ftrs),
      mz_pattern,
      mz_replacement )

    feature_cols <- c("mz", "time", "mz_min", "mz_max")
    sample_cols_idx <- which(!colnames(final.ftrs) %in% feature_cols)

    colnames(final.ftrs)[sample_cols_idx] <- tools::file_path_sans_ext(colnames(final.ftrs)[sample_cols_idx])
    colnames(final.times)[sample_cols_idx] <- tools::file_path_sans_ext(colnames(final.times)[sample_cols_idx])

    sample_cols <- colnames(final.ftrs)[sample_cols_idx]

    final.ftrs <- tidyr::pivot_longer(final.ftrs,
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "sample_intensity"
    )

    final.times <- tidyr::pivot_longer(final.times,
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "sample_rt"
    )

    recovered_feature_sample_table <- dplyr::full_join(final.times, final.ftrs,
      by = c("mz", "time", "mz_min", "mz_max", "sample")
    )

    colnames(recovered_feature_sample_table)[2] <- "rt"

    batchwise[[batch]]$recovered_feature_sample_table <- recovered_feature_sample_table

    batchwise[[batch]]$corrected_features <- batchwise[[batch]]$features2
    batchwise[[batch]]$extracted_features <- batchwise[[batch]]$features
  }

  return(batchwise)
}

two.step.hybrid <- function(filenames,
                            metadata,
                            work_dir,
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
                            BIC.factor = 2) {
  filenames_batchwise <- bind_batch_label_column(filenames, metadata)
  batches_idx <- unique(metadata$batch)
  batchwise <- new("list")
  message("* processing ", length(batches_idx), " batches separately")

  for (batch.i in batches_idx) {
    files_batch <- dplyr::filter(filenames_batchwise, batch == batch.i)$filename
    message("* processing ", length(files_batch), " samples from batch ", batch.i)

    features <- semi.sup(
      files = files_batch,
      folder = work_dir,
      n.nodes = cluster,
      known.table = known.table,
      sd.cut = sd.cut,
      sigma.ratio.lim = sigma.ratio.lim,
      component.eliminate = component.eliminate,
      moment.power = moment.power,
      min.pres = min.pres,
      min.run = min.run,
      min.exp = ceiling(min.within.batch.prop.detect * length(files_batch)),
      mz.tol = mz.tol,
      baseline.correct.noise.percentile = baseline.correct.noise.percentile,
      baseline.correct = baseline.correct,
      align.mz.tol = align.mz.tol,
      align.chr.tol = align.chr.tol,
      max.align.mz.diff = max.align.mz.diff,
      recover.mz.range = recover.mz.range,
      recover.chr.range = recover.chr.range,
      use.observed.range = use.observed.range,
      shape.model = shape.model,
      new.feature.min.count = new.feature.min.count,
      recover.min.count = recover.min.count
    )

    features$final.ftrs <- features$final.ftrs[order(features$final.ftrs[, 1], features$final.ftrs[, 2]), ]
    features$final.times <- features$final.times[order(features$final.times[, 1], features$final.times[, 2]), ]

    batchwise[[batch.i]] <- features
  }

  batchwise <- semisup_to_hybrid_adapter(batchwise, batches_idx)

  step_one_features <- list()
  for (batch_id in batches_idx) {
    step_one_features[[batch_id]] <- compute_intensity_medians(
      batchwise[[batch_id]]$recovered_feature_sample_table
    ) %>%
      dplyr::select(-c("mz_min", "mz_max"))
  }

  cluster <- parallel::makeCluster(cluster)
  doParallel::registerDoParallel(cluster)
  clusterEvalQ(cluster, library(recetox.aplcms))

  extracted_features <- list()
  for (batch_id in batches_idx) {
    extracted <- batchwise[[batch_id]]$recovered_feature_sample_table
    extracted <- long_to_wide_feature_table(extracted)
    intensities <- dplyr::select(extracted, contains("_intensity"))
    median_intensity <- apply(intensities, 1, median)
    extracted <- dplyr::select(extracted, mz, rt)
    extracted[, 3:4] <- NA
    colnames(extracted)[3:4] <- c("mz_min", "mz_max")
    extracted[, 5] <- median_intensity
    colnames(extracted)[5] <- "median_intensity"
    extracted_features[[batch_id]] <- extracted
  }

  message("* aligning time")
  corrected <- adjust.time(extracted_features,
    mz.tol = batch.align.mz.tol,
    chr.tol = batch.align.chr.tol,
    find.tol.max.d = 10 * mz.tol,
    max.align.mz.diff = max.align.mz.diff,
    rt_colname = "rt"
  )

  message("* aligning features")
  aligned <- align_features(
    sample_names = paste0("batch_", batches_idx),
    features = corrected,
    min.exp = ceiling(min.batch.prop * length(batches_idx)),
    mz.tol = batch.align.mz.tol,
    chr.tol = batch.align.chr.tol,
    find.tol.max.d = 10 * mz.tol,
    max.align.mz.diff = max.align.mz.diff,
    rt_colname = "rt"
  )

  aligned_wide <- as_wide_aligned_table(aligned)

  message("* recovering features across batches")
  recovered <- feature_recovery(
    cluster = cluster,
    step_one_features = step_one_features,
    batchwise = batchwise,
    filenames_batchwise = filenames_batchwise,
    corrected = corrected,
    aligned = aligned_wide,
    batches_idx = batches_idx,
    mz.tol = mz.tol,
    batch.align.mz.tol = batch.align.mz.tol,
    batch.align.chr.tol = batch.align.chr.tol,
    recover.mz.range = recover.mz.range,
    recover.chr.range = recover.chr.range,
    use.observed.range = use.observed.range,
    min.bw = min.bw,
    max.bw = max.bw,
    recover.min.count = recover.min.count
  )

  stopCluster(cluster)

  recovered_filtered <- filter_features_by_presence(
    feature_table = recovered,
    metadata = metadata,
    batches_idx = batches_idx,
    within_batch_threshold = min.within.batch.prop.report,
    across_batch_threshold = min.batch.prop
  )

  recovered_filtered <- dplyr::arrange(recovered_filtered, mz, rt)
  recovered_features <- wide_to_long_feature_table(recovered_filtered)

  features <- new("list")
  features$batchwise_features <- batchwise
  features$aligned_features <- as_feature_sample_table(
    rt_crosstab = aligned$rt_crosstab,
    int_crosstab = aligned$int_crosstab
  )
  features$corrected_features <- corrected
  features$final_features <- recovered_features
  return(features)
}
