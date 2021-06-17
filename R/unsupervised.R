as_feature_crosstab <- function(feature_names, sample_names, data) {
  colnames(data) <- c('mz', 'rt', 'mz_min', 'mz_max', sample_names)
  rownames(data) <- feature_names
  as.data.frame(data)
}

as_feature_sample_table <- function(rt_crosstab, int_crosstab) {
  feature_names <- rownames(rt_crosstab)
  sample_names <- colnames(rt_crosstab)[-(1:4)]

  feature_table <- data.frame(
    feature = feature_names,
    mz = rt_crosstab[, 1],
    rt = rt_crosstab[, 2]
  )

  # series of conversions to produce a table type from data.frame
  rt_crosstab <- as.table(as.matrix(rt_crosstab[, -(1:4)]))
  int_crosstab <- as.table(as.matrix(int_crosstab[, -(1:4)]))

  crosstab_axes <- list(feature = feature_names, sample = sample_names)
  dimnames(rt_crosstab) <- dimnames(int_crosstab) <- crosstab_axes

  x <- as.data.frame(rt_crosstab, responseName = 'sample_rt')
  y <- as.data.frame(int_crosstab, responseName = 'sample_intensity')

  data <- merge(x, y, by = c('feature', 'sample'))
  data <- merge(feature_table, data, by = 'feature')
  data
}

check_files <- function(filenames) {
  missing <- !file.exists(filenames)
  missing_filenames <- paste0('\t', filenames[missing], collapse = '\n')

  if (any(missing)) {
    stop("Cannot find the following files:\n", missing_filenames)
  }
}

get_sample_name <- function(filename) {
  tools::file_path_sans_ext(basename(filename))
}

sort_samples_by_acquisition_number <- function (filenames) {
  # assumes that the filenames contain an acquisition number
  # ideal solution would be to read the acquisition number directly from mzml
  sort(unlist(filenames))
}

extract_features <- function(
  cluster,
  filenames,
  min_pres,
  min_run,
  mz_tol,
  baseline_correct,
  baseline_correct_noise_percentile,
  intensity_weighted,
  min_bandwidth,
  max_bandwidth,
  sd_cut,
  sigma_ratio_lim,
  shape_model,
  peak_estim_method,
  component_eliminate,
  moment_power,
  BIC_factor
) {
  clusterExport(cluster, list(
    'proc.cdf',
    'prof.to.features',
    'load.lcms',
    'adaptive.bin',
    'find.turn.point',
    'merge.seq.3',
    'cont.index',
    'interpol.area'
  ))

  parLapply(cluster, filenames, function(filename) {
    profile <- proc.cdf(
      filename = filename,
      min.pres = min_pres,
      min.run = min_run,
      tol = mz_tol,
      baseline.correct = baseline_correct,
      baseline.correct.noise.percentile = baseline_correct_noise_percentile,
      intensity.weighted = intensity_weighted,
      do.plot = FALSE,
      cache = FALSE
    )
    features <- prof.to.features(
      a = profile,
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
}

align_features <- function(sample_names, ...) {
  aligned <- feature.align(...)
  feature_names <- seq_len(nrow(aligned$pk.times))

  list(
    mz_tolerance = as.numeric(aligned$mz.tol),
    rt_tolerance = as.numeric(aligned$chr.tol),
    rt_crosstab = as_feature_crosstab(feature_names, sample_names, aligned$pk.times),
    int_crosstab = as_feature_crosstab(feature_names, sample_names, aligned$aligned.ftrs)
  )
}

recover_weaker_signals <- function(
  cluster,
  filenames,
  extracted_features,
  corrected_features,
  aligned_rt_crosstab,
  aligned_int_crosstab,
  original_mz_tolerance,
  aligned_mz_tolerance,
  aligned_rt_tolerance,
  mz_range,
  rt_range,
  use_observed_range,
  min_bandwidth,
  max_bandwidth,
  recover_min_count
) {
  clusterExport(cluster, 'recover.weaker')

  recovered <- parLapply(cluster, seq_along(filenames), function(i) {
    recover.weaker(
      loc = i,
      filename = filenames[[i]],
      this.f1 = extracted_features[[i]],
      this.f2 = corrected_features[[i]],
      pk.times = aligned_rt_crosstab,
      aligned.ftrs = aligned_int_crosstab,
      orig.tol = original_mz_tolerance,
      align.mz.tol = aligned_mz_tolerance,
      align.chr.tol = aligned_rt_tolerance,
      mz.range = mz_range,
      chr.range = rt_range,
      use.observed.range = use_observed_range,
      bandwidth = 0.5,
      min.bw = min_bandwidth,
      max.bw = max_bandwidth,
      recover.min.count = recover_min_count
    )
  })

  feature_table <- aligned_rt_crosstab[, 1:4]
  rt_crosstab <- cbind(feature_table, sapply(recovered, function(x) x$this.times))
  int_crosstab <- cbind(feature_table, sapply(recovered, function(x) x$this.ftrs))

  feature_names <- rownames(feature_table)
  sample_names <- colnames(aligned_rt_crosstab[, -(1:4)])

  list(
    extracted_features = lapply(recovered, function(x) x$this.f1),
    corrected_features = lapply(recovered, function(x) x$this.f2),
    rt_crosstab = as_feature_crosstab(feature_names, sample_names, rt_crosstab),
    int_crosstab = as_feature_crosstab(feature_names, sample_names, int_crosstab)
  )
}

unsupervised <- function(
  filenames,
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
  recover_mz_range = NA,
  recover_chr_range = NA,
  use_observed_range = TRUE,
  recover_min_count = 3,
  intensity_weighted = FALSE,
  cluster = parallel::detectCores()
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

  message("**** weaker signal recovery ****")
  recovered <- recover_weaker_signals(
    cluster = cluster,
    filenames = filenames,
    extracted_features = extracted,
    corrected_features = corrected,
    aligned_rt_crosstab = aligned$rt_crosstab,
    aligned_int_crosstab = aligned$int_crosstab,
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

  aligned_feature_sample_table <- as_feature_sample_table(
    rt_crosstab = aligned$rt_crosstab,
    int_crosstab = aligned$int_crosstab
  )
  recovered_feature_sample_table <- as_feature_sample_table(
    rt_crosstab = recovered$rt_crosstab,
    int_crosstab = recovered$int_crosstab
  )

  list(
    extracted_features = recovered$extracted_features,
    corrected_features = recovered$corrected_features,
    aligned_feature_sample_table = aligned_feature_sample_table,
    recovered_feature_sample_table = recovered_feature_sample_table,
    aligned_mz_toletance = as.numeric(aligned$mz_tolerance),
    aligned_rt_tolerance = as.numeric(aligned$rt_tolerance)
  )
}
