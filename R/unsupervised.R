.check_tasks <- function(message, tasks, results) {
  errors <- vapply(results, function (x) inherits(x, "try-error"), logical(1))
  
  if (any(errors)) {
    stop(
      message,
      ":\n",
      paste("\t", tasks[errors], ":", results[errors], collapse = "\n")
    )
  }
  
  return(results)
}

.check_files <- function(files) {
  if (any(missing <- !file.exists(files))) {
    stop(
      "Cannot find the following files:\n",
      paste0("\t", files[missing], collapse = "\n")
    )
  }
}

extract_features <- function(
  files,
  min_pres,
  min_run,
  mz_tol,
  baseline_correct_noise_percentile,
  shape_model,
  BIC_factor,
  baseline_correct,
  peak_estim_method,
  min_bandwidth,
  max_bandwidth,
  sd_cut,
  sigma_ratio_lim,
  component_eliminate,
  moment_power,
  intensity_weighted,
  cluster = NULL
) {
  .check_files(files)

  if (is.null(cluster)) {
    cluster <- parallel::makeCluster(parallel::detectCores())
    on.exit(parallel::stopCluster(cluster))
  }

  features <- parallel::parLapply(cluster, files, function(file) {
    try({
      prof <- proc.cdf(
        filename = file,
        min.pres = min_pres,
        min.run = min_run,
        tol = mz_tol,
        baseline.correct = baseline_correct,
        baseline.correct.noise.percentile = baseline_correct_noise_percentile,
        do.plot = FALSE,
        intensity.weighted = intensity_weighted,
        cache = FALSE
      )
      feat <- prof.to.features(
        a = prof,
        min.bw = min_bandwidth,
        max.bw = max_bandwidth,
        sd.cut = sd_cut,
        sigma.ratio.lim = sigma_ratio_lim,
        shape.model = shape_model,
        estim.method = peak_estim_method,
        do.plot = FALSE,
        component.eliminate = component_eliminate,
        power = moment_power,
        BIC.factor = BIC_factor
      )
    })
  })

  .check_tasks("Feature extraction was unsuccesfull", files, features)
}

recover_weaker_signals <- function(
  files,
  features,
  corrected_features,
  aligned_features,
  pk_times,
  aligned_mz_tol,
  aligned_chr_tol,
  mz_range,
  chr_range,
  use_observed_range,
  orig_tol,
  min_bandwidth,
  max_bandwidth,
  recover_min_count,
  cluster = NULL
) {
  .check_files(files)

  if (is.null(cluster)) {
    cluster <- parallel::makeCluster(parallel::detectCores())
    on.exit(parallel::stopCluster(cluster))
  }

  recovered <- parallel::parLapply(cluster, seq_along(files), function(i) {
    try({
      recover.weaker(
        filename = files[[i]],
        loc = i,
        aligned.ftrs = aligned_features,
        pk.times = pk_times,
        align.mz.tol = aligned_mz_tol,
        align.chr.tol = aligned_chr_tol,
        this.f1 = features[[i]],
        this.f2 = corrected_features[[i]],
        mz.range = mz_range,
        chr.range = chr_range,
        use.observed.range = use_observed_range,
        orig_tol,
        min.bw = min_bandwidth,
        max.bw = max_bandwidth,
        bandwidth = 0.5,
        recover.min.count = recover_min_count
      )
    })
  })

  .check_tasks("Signal recovery was unsuccesfull", files, recovered)

  list(
    f1 = lapply(recovered, function(x) x$this.f1),
    f2 = lapply(recovered, function(x) x$this.f2),
    times = simplify2array(lapply(recovered, function(x) x$this.times)),
    features = simplify2array(lapply(recovered, function(x) x$this.ftrs))
  )
}

unsupervised <- function(
  files,
  min_exp = 2,
  min_pres = 0.5,
  min_run = 12,
  mz_tol = 1e-05,
  baseline_correct_noise_percentile = 0.05,
  shape_model = "bi-Gaussian",
  BIC_factor = 2,
  baseline_correct = 0,
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
  recover_mz_range = NA,
  recover_chr_range = NA,
  use_observed_range = TRUE,
  recover_min_count = 3,
  intensity_weighted = FALSE,
  cluster = parallel::detectCores()
) {
  if (is.numeric(cluster)) {
    cluster <- parallel::makeCluster(cluster)
    on.exit(parallel::stopCluster(cluster))
  } else if (!is(cluster, "cluster")) {
    stop("unsupported value for `cluster` parameter: ", cluster)
  }

  # NOTE: side effect (doParallel has no functionality to clean up)
  doParallel::registerDoParallel(cluster)

  # further processing requires sorted file list
  files <- sort(unlist(files))

  message("**** feature extraction ****")
  features <- extract_features(
    files = files,
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
    BIC_factor = BIC_factor,
    cluster = cluster
  )

  message("**** time correction ****")
  corrected <- adjust.time(
    features = features,
    mz.tol = align_mz_tol,
    chr.tol = align_chr_tol,
    find.tol.max.d = 10 * mz_tol,
    max.align.mz.diff = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** feature alignemnt ****")
  aligned <- feature.align(
    features = corrected,
    min.exp = min_exp,
    mz.tol = align_mz_tol,
    chr.tol = align_chr_tol,
    find.tol.max.d = 10 * mz_tol,
    max.align.mz.diff = max_align_mz_diff,
    do.plot = FALSE
  )

  message("**** weaker signal recovery ****")
  recovered <- recover_weaker_signals(
    files = files,
    features = features,
    corrected_features = corrected,
    aligned_features = aligned$aligned.ftrs,
    pk_times = aligned$pk.times,
    aligned_mz_tol = aligned$mz.tol,
    aligned_chr_tol = aligned$chr.tol,
    mz_range = recover_mz_range,
    chr_range = recover_chr_range,
    use_observed_range = use_observed_range,
    orig_tol = mz_tol,
    min_bandwidth = min_bandwidth,
    max_bandwidth = max_bandwidth,
    recover_min_count = recover_min_count,
    cluster = cluster
  )

  recovered_times <- cbind(aligned$pk.times[, 1:4], recovered$times)
  recovered_features <- cbind(aligned$aligned.ftrs[, 1:4], recovered$features)

  colnames(recovered_times) <-
    colnames(aligned$pk.times) <-
      c("mz", "rt", "mz_min", "mz_max", paste0("time.", basename(files)))
  colnames(aligned$aligned.ftrs) <-
    colnames(recovered_features) <-
      c("mz", "rt", "mz_min", "mz_max", paste0("intensity.", basename(files)))

  aligned_times <- subset(aligned$pk.times, select = -c(mz_min, mz_max))
  recovered_times <- subset(recovered_times, select = -c(mz_min, mz_max))

  final_peaks <- merge(recovered_features, recovered_times, by = c("mz", "rt"))
  aligned_peaks <- merge(aligned$aligned.ftrs, aligned_times, by = c("mz", "rt"))

  list(
    final_peaks = as.data.frame(final_peaks),
    aligned_peaks = as.data.frame(aligned_peaks),
    aligned_mz_tolerance = as.numeric(aligned$mz.tol),
    aligned_rt_tolerance = as.numeric(aligned$chr.tol),
    extracted_features = as.data.frame(do.call(rbind, recovered$f1)),
    corrected_features = as.data.frame(do.call(rbind, recovered$f2))
  )
}
