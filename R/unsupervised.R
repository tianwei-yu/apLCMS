#' @import tools snow splines parallel doParallel
NULL
#> NULL

#' @importFrom dplyr select inner_join
as_feature_crosstab <- function(sample_names, metadata, data) {
  metadata_cols <- c('id', 'mz', 'rt', 'mzmin', 'mzmax')
  data <- select(metadata, metadata_cols) |>
    inner_join(data, on='id')
  colnames(data) <- c(metadata_cols, sample_names)

  return(data)
}

as_feature_sample_table <- function(metadata, rt_crosstab, int_crosstab) {
  feature_names <- as.character(rt_crosstab$id)
  sample_names <- colnames(metadata)[-c(1:8)]

  feature_table <- data.frame(
    feature = feature_names,
    mz = metadata$mz,
    rt = metadata$rt
  )

  # series of conversions to produce a table type from data.frame
  rt_crosstab <- as.table(as.matrix(rt_crosstab[, -1]))
  int_crosstab <- as.table(as.matrix(int_crosstab[, -1]))

  crosstab_axes <- list(feature = feature_names, sample = sample_names)
  dimnames(rt_crosstab) <- dimnames(int_crosstab) <- crosstab_axes

  x <- as.data.frame(rt_crosstab, responseName = 'sample_rt')
  y <- as.data.frame(int_crosstab, responseName = 'sample_intensity')

  data <- merge(x, y, by = c('feature', 'sample'))
  data <- merge(feature_table, data, by = 'feature')
  return(data)
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
  recover_mz_range,
  recover_rt_range,
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
      recover_mz_range = recover_mz_range,
      recover_rt_range = recover_rt_range,
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
#' @param baseline_correct_noise_percentile The percentile of signal strength of those EIC that don't pass the 
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
#' @param mz_tol_relative The m/z tolerance level for peak alignment. The default is NA, which allows the program to search 
#'  for the tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, 
#'  multiplied by the m/z value, becomes the cutoff level.
#' @param rt_tol_relative The retention time tolerance level for peak alignment. The default is NA, which allows the program 
#'  to search for the tolerance level based on the data.
#' @param mz_tol_absolute As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable when the 
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
  bandwidth = 0.5,
  min_bandwidth = NA,
  max_bandwidth = NA,
  sd_cut = c(0.01, 500),
  sigma_ratio_lim = c(0.01, 100),
  component_eliminate = 0.01,
  moment_power = 1,
  mz_tol_relative = NA,
  rt_tol_relative = NA,
  mz_tol_absolute = 0.01,
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
          min_pres = min_pres,
          min_run = min_run,
          mz_tol = mz_tol,
          baseline_correct = baseline_correct,
          baseline_correct_noise_percentile = baseline_correct_noise_percentile,
          intensity_weighted = intensity_weighted,
          do.plot = FALSE,
          cache = FALSE
      )
  })
  
  feature_tables <- snow::parLapply(cluster, profiles, function(profile) {
      prof.to.features(
          profile = profile,
          bandwidth = bandwidth,
          min_bandwidth = min_bandwidth,
          max_bandwidth = max_bandwidth,
          sd_cut = sd_cut,
          sigma_ratio_lim = sigma_ratio_lim,
          shape_model = shape_model,
          peak_estim_method = peak_estim_method,
          component_eliminate = component_eliminate,
          moment_power = moment_power,
          BIC_factor = BIC_factor,
          do.plot = FALSE
      )
  })

  message("**** computing clusters ****")
  extracted_clusters <- compute_clusters(
    feature_tables = feature_tables,
    mz_tol_relative = mz_tol_relative,
    mz_tol_absolute = mz_tol_absolute,
    mz_max_diff = 10 * mz_tol,
    rt_tol_relative = rt_tol_relative,
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
    rt_tol_relative = rt_tol_relative
  )

  message("**** feature alignment ****")
  aligned <- create_aligned_feature_table(
      dplyr::bind_rows(adjusted_clusters$feature_tables),
      min_exp,
      sample_names,
      adjusted_clusters$rt_tol_relative,
      adjusted_clusters$mz_tol_relative
  )

  message("**** weaker signal recovery ****")
  recovered <- lapply(seq_along(filenames), function(i) {
    recover.weaker(
      filename = filenames[[i]],
      sample_name = as.character(i),
      extracted_features = feature_tables[[i]],
      adjusted_features = corrected[[i]],
      metadata_table = aligned$metadata,
      rt_table = aligned$rt,
      intensity_table = aligned$intensity,
      orig.tol = mz_tol,
      align.mz.tol = aligned$mz_tol_relative,
      align.rt.tol = aligned$rt_tol_relative,
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

  message("**** computing clusters ****")
  recovered_clusters <- compute_clusters(
    feature_tables = recovered_adjusted,
    mz_tol_relative = adjusted_clusters$mz_tol_relative,
    mz_tol_absolute = adjusted_clusters$rt_tol_relative,
    mz_max_diff = 10 * mz_tol,
    rt_tol_relative = rt_tol_relative
  )

  message("**** feature alignment ****")
  recovered_aligned <- create_aligned_feature_table(
      dplyr::bind_rows(recovered_clusters$feature_tables),
      min_exp,
      sample_names,
      recovered_clusters$rt_tol_relative,
      recovered_clusters$mz_tol_relative
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
    corrected_features = recovered$adjusted_features,
    aligned_feature_sample_table = aligned_feature_sample_table,
    recovered_feature_sample_table = recovered_feature_sample_table,
    aligned_mz_tolerance = as.numeric(aligned$mz_tol_relative),
    aligned_rt_tolerance = as.numeric(aligned$rt_tol_relative)
  )
}
