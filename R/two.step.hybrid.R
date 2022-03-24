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
    dplyr::inner_join(sample_rts, on = c("mz", "rt")) %>%
    dplyr::inner_join(sample_intensities, on = c("mz", "rt"))
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
  filenames <- inner_join(filenames, metadata, on = "sample_name")
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

  aligned <- as_feature_sample_table(
    rt_crosstab = aligned$rt_crosstab,
    int_crosstab = aligned$int_crosstab
  )

  stopCluster(cl)

  message("*** recovering features across batches ***")

  cl <- makeCluster(cluster)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(recetox.aplcms))

  for (batch_id in batches_idx)
  {
    this.fake <- step_one_features[[batch_id]] # intensity per sample table after step one
    this.fake.time <- batchwise[[batch_id]]$final.times # rt per sample table after step one
    # this.medians <- apply(this.fake[, -1:-4], 1, median) # intensity median, commented out because already present in this.fake

    this.pk.time <- this.aligned <- matrix(0, nrow = nrow(aligned), ncol = ncol(this.fake) - 4) # zero matrix with dimensions (num_features x num_samples)
    # adjusting the time (already within batch adjusted)
    this.features <- readjust_times(batchwise[[batch_id]], corrected[[batch_id]])

    for (i in 1:nrow(this.aligned))
    {
      if (fake3$aligned[i, batch_id + 4] != 0) {
        sel <- which(fake3$aligned[i, 3] <= this.fake[, 1] & fake3$aligned[i, 4] >= this.fake[, 1] & abs(this.medians - fake3$aligned[i, batch_id + 4]) < 1)
        if (length(sel) == 0) sel <- which(fake3$aligned[i, 3] <= this.fake[, 1] & fake3$aligned[i, 4] >= this.fake[, 1])
        if (length(sel) == 0) {
          message("batch", batch_id, " row ", i, " match issue.")
        } else {
          if (length(sel) == 1) {
            this.aligned[i, ] <- this.fake[sel, -1:-4]
            this.pk.time[i, ] <- this.fake.time[sel, -1:-4]
          } else {
            this.aligned[i, ] <- apply(this.fake[sel, -1:-4], 2, sum)
            this.pk.time[i, ] <- apply(this.fake.time[sel, -1:-4], 2, median)
          }
        }
      } else {
        ### go into individual feature tables to find a match
        recaptured <- rep(0, ncol(this.aligned))
        recaptured.time <- rep(NA, ncol(this.aligned))
        for (j in 1:length(this.features))
        {
          diff.mz <- abs(this.features[[j]][, 1] - fake3$aligned[i, 1])
          diff.time <- abs(this.features[[j]][, 2] - fake3$aligned[i, 2])
          sel <- which(diff.mz < fake3$aligned[i, 1] * batch.align.mz.tol & diff.time <= batch.align.chr.tol)
          if (length(sel) > 1) sel <- sel[which(diff.time[sel] == min(diff.time[sel]))[1]]

          if (length(sel) == 1) {
            recaptured[j] <- this.features[[j]][sel, 5]
            recaptured.time[j] <- this.features[[j]][sel, 2]
          }
        }
        this.aligned[i, ] <- recaptured
        this.pk.time[i, ] <- recaptured.time
      }
    }

    colnames(this.aligned) <- colnames(this.fake)[-1:-4]
    colnames(this.pk.time) <- colnames(this.fake)[-1:-4]

    #### go back to cdf files

    that.aligned <- cbind(fake3$aligned.ftrs[, 1:4], this.aligned)
    that.pk.time <- cbind(fake3$aligned.ftrs[, 1:4], this.pk.time)

    # for(i in 1:ncol(this.aligned))
    new.this.aligned <- foreach(i = 1:ncol(this.aligned), .combine = cbind) %dopar% {
      r <- recover.weaker(filename = colnames(this.aligned)[i], loc = i, aligned.ftrs = that.aligned, pk.times = that.pk.time, align.mz.tol = batch.align.mz.tol, align.chr.tol = batch.align.chr.tol, this.f1 = batchwise[[batch_id]]$features[[i]], this.f2 = batchwise[[batch_id]]$features2[[i]], mz.range = recover.mz.range, chr.range = recover.chr.range, use.observed.range = use.observed.range, orig.tol = mz.tol, min.bw = min.bw, max.bw = max.bw, bandwidth = .5, recover.min.count = recover.min.count)

      r$this.ftrs
    }

    colnames(new.this.aligned) <- colnames(this.aligned)

    if (batch_id == 1) {
      aligned <- new.this.aligned
    } else {
      aligned <- cbind(aligned, new.this.aligned)
    }
  }

  aligned <- cbind(fake3$aligned.ftrs[, 1:4], aligned)

  batch.presence.mat <- matrix(0, nrow = nrow(aligned), ncol = length(batches))
  for (batch_id in 1:length(batches))
  {
    batch <- batches[batch_id]
    this.mat <- aligned[, which(colnames(aligned) %in% info[info[, 2] == batch, 1])]
    this.mat <- 1 * (this.mat != 0)
    this.presence <- apply(this.mat, 1, sum) / ncol(this.mat)
    batch.presence.mat[, batch_id] <- 1 * (this.presence >= min.within.batch.prop.report)
  }
  batch.presence <- apply(batch.presence.mat, 1, sum) / ncol(batch.presence.mat)
  final.aligned <- aligned[which(batch.presence >= min.batch.prop), ]

  stopCluster(cl)
  to.return <- new("list")
  to.return$batchwise.results <- batchwise
  to.return$all.detected.ftrs <- aligned
  to.return$final.ftrs <- final.aligned
  return(to.return)
}
