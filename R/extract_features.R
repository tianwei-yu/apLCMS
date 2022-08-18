#' @import snow
NULL
#> NULL

#' feature extraction
#' 
#' extract features
#' 
#' @param cluster The number of CPU cores to be used
#' @param filenames The CDF file names.
#' @param min_pres This is a parameter of thr run filter, to be passed to the function proc.cdf().
#' @param min_run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z 
#'  to be considered a peak.
#' @param mz_tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of 
#'  the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is 
#'  the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline_correct After grouping the observations, the highest intensity in each group is found. If the highest is lower than 
#'  this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the 
#'  height of the noise groups.
#' @param baseline_correct_noise_percentile The perenctile of signal strength of those EIC that don't pass the 
#'  run filter, to be used as the baseline threshold of signal strength.
#' @param intensity_weighted Whether to use intensity to weight mass density estimation.
#' @param min_bandwidth The minimum bandwidth to use in the kernel smoother.
#' @param max_bandwidth The maximum bandwidth to use in the kernel smoother.
#' @param sd_cut A vector of two. Features with standard deviation outside the range defined by the two numbers 
#'  are eliminated.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard 
#'  deviation and the righ-standard deviation of the bi-Gaussian fuction used to fit the data.
#' @param shape_model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and 
#'  "Gaussian". When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param peak_estim_method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param component_eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts 
#'  for a proportion of intensities less than this value, the component will be ignored.
#' @param moment_power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture 
#'  model in an EIC.
#' @param BIC_factor The factor that is multiplied on the number of parameters to modify the BIC criterion. 
#'  If larger than 1, models with more peaks are penalized more.
#' @export
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
  snow::clusterExport(cluster, list(
    'proc.cdf',
    'prof.to.features',
    'load.lcms',
    'adaptive.bin',
    'find.turn.point',
    'msExtrema',
    'find_local_maxima',
    'combine.seq.3',
    'cont.index',
    'interpol.area',
    'load_file',
    'load_data',
    'plot_raw_profile_histogram',
    'compute_mass_values',
    'compute_densities',
    'compute_breaks',
    'compute_breaks_3',
    'compute_boundaries',
    'increment_counter',
    'rm.ridge',
    'compute_delta_rt',
    'bigauss.mix',
    'bigauss.esti',
    'rev_cum_sum',
    'compute_bounds',
    'validate_inputs',
    'preprocess_bandwidth',
    'preprocess_profile',
    'compute_gaussian_peak_shape',
    'compute_chromatographic_profile',
    'compute_dx',
    'compute_initiation_params',
    'compute_e_step',
    'compute_start_bound',
    'compute_end_bound',
    'compute_bounds',
    'compute_scale'
  ))
  snow::clusterEvalQ(cluster, library("dplyr"))


  snow::parLapply(cluster, filenames, function(filename) {
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
}
