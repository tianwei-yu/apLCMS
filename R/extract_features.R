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
    'combine.seq.3',
    'cont.index',
    'interpol.area',
    'load_file',
    'load_data',
    'plot_raw_profile_histogram'
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
