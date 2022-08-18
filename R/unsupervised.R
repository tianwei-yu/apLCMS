#' @import tools snow splines parallel doParallel
NULL
#> NULL

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

#' @export
sort_samples_by_acquisition_number <- function (filenames) {
  # assumes that the filenames contain an acquisition number
  # ideal solution would be to read the acquisition number directly from mzml
  sort(unlist(filenames))
}

align_features <- function(sample_names, ...) {
  aligned <- feature.align(...)
  feature_names <- seq_len(nrow(aligned$peak_times))

  list(
    mz_tolerance = as.numeric(aligned$mz_tol_relative),
    rt_tolerance = as.numeric(aligned$rt_tol_relative),
    rt_crosstab = as_feature_crosstab(feature_names, sample_names, aligned$peak_times),
    int_crosstab = as_feature_crosstab(feature_names, sample_names, aligned$aligned_features)
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
  snow::clusterExport(cluster, c('recover.weaker'))
  snow::clusterEvalQ(cluster, library("splines"))

  recovered <- lapply(seq_along(filenames), function(i) {
    recover.weaker(
      sample_name = get_sample_name(filenames[i]),
      filename = filenames[[i]],
      extracted_features = as_tibble(extracted_features[[i]]),
      adjusted_features = as_tibble(corrected_features[[i]]),
      pk.times = aligned_rt_crosstab,
      aligned.ftrs = aligned_int_crosstab,
      orig.tol = original_mz_tolerance,
      align.mz.tol = aligned_mz_tolerance,
      align.rt.tol = aligned_rt_tolerance,
      mz.range = mz_range,
      rt.range = rt_range,
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

#' Runs features extraction in unsupervised mode.
#' 
#' features extraction in unsupervised mode.
#' 
#' @param filenames The CDF file names.
#' @param min_exp A feature has to show up in at least this number of profiles to be included in the final result.
#' @param min_pres This is a parameter of the run filter, to be passed to the function proc.cdf().
#' @param min_run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z 
#'  to be considered a peak.
#' @param mz_tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of 
#'  the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is 
#'  the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline_correct After grouping the observations, the highest intensity in each group is found. 
#'  If the highest is lower than this value, the entire group will be deleted. The default value is NA, in 
#'  which case the program uses the 75th percentile of the height of the noise groups.
#' @param baseline_correct_noise_percentile The perenctile of signal strength of those EIC that don't pass the 
#'  run filter, to be used as the baseline threshold of signal strength.
#' @param shape_model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and 
#'  "Gaussian". When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param BIC_factor The factor that is multiplied on the number of parameters to modify the BIC criterion. 
#'  If larger than 1, models with more peaks are penalized more.
#' @param peak_estim_method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param min_bandwidth The minimum bandwidth to use in the kernel smoother.
#' @param max_bandwidth The maximum bandwidth to use in the kernel smoother.
#' @param sd_cut A vector of two. Features with standard deviation outside the range defined by the two numbers 
#'  are eliminated.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard 
#'  deviation and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @param component_eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts 
#'  for a proportion of intensities less than this value, the component will be ignored.
#' @param moment_power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture 
#'  model in an EIC.
#' @param align_mz_tol The m/z tolerance level for peak alignment. The default is NA, which allows the program to search 
#'  for the tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, 
#'  multiplied by the m/z value, becomes the cutoff level.
#' @param align_rt_tol The retention time tolerance level for peak alignment. The default is NA, which allows the program 
#'  to search for the tolerance level based on the data.
#' @param max_align_mz_diff As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable when the 
#'  m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly influences feature matching 
#'  in higher m/z range.
#' @param recover_mz_range The m/z around the feature m/z to search for observations. The default value is NA, in which 
#'  case 1.5 times the m/z tolerance in the aligned object will be used.
#' @param recover_rt_range The retention time around the feature retention time to search for observations. The default 
#'  value is NA, in which case 0.5 times the retention time tolerance in the aligned object will be used.
#' @param use_observed_range If the value is TRUE, the actual range of the observed locations of the feature in all 
#'  the spectra will be used.
#' @param recover_min_count Minimum number of raw data points to support a recovery.
#' @param intensity_weighted Whether to use intensity to weight mass density estimation.
#' @param cluster The number of CPU cores to be used
#' @export
#' @examples
#' unsupervised(test_files, cluster = num_workers)
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
  align_rt_tol = NA,
  max_align_mz_diff = 0.01,
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

  message("**** feature alignemnt ****")
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
    rt_range = recover_rt_range,
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
    aligned_mz_tolerance = as.numeric(aligned$mz_tolerance),
    aligned_rt_tolerance = as.numeric(aligned$rt_tolerance)
  )
}
