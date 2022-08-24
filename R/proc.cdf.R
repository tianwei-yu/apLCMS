#' @import tidyr
NULL
#> NULL

#' Load raw data from file
#' @export
load_file <- function(filename) {
  this <- load.lcms(filename)
  this <- tidyr::drop_na(this)
  return(this)
}

#' Load data either from cache or load raw file and detect peaks.
#' @export
load_data <- function(filename,
                      cache,
                      min_elution_length,
                      min_presence,
                      mz_tol,
                      baseline_correct,
                      intensity_weighted) {
  rawprof_filename <- paste(strsplit(tolower(filename), "\\.")[[1]][1], "_", min_elution_length, "_", min_presence, "_", mz_tol, ".rawprof", sep = "")

  if (cache && file.exists(rawprof_filename)) {
    load(rawprof_filename)
  } else {
    raw.data <- load_file(filename)
    raw.prof <- adaptive.bin(
      raw.data,
      min_elution_length = min_elution_length,
      min_presence = min_presence,
      mz_tol = mz_tol,
      baseline_correct = baseline_correct,
      intensity_weighted = intensity_weighted
    )
  }

  if (cache && !file.exists(rawprof_filename)) {
    save(raw.prof, file = rawprof_filename)
  }

  return(raw.prof)
}

#' Filter noise and detect peaks from LC/MS data in CDF format
#' 
#' This function applies the run filter to remove noise. Data points are grouped into EICs in this step.
#' 
#' @param filename The CDF file name. If the file is not in the working directory, the path needs to be given.
#' @param min_presence Run filter parameter. The minimum proportion of presence in the time period for a series of 
#'  signals grouped by m/z to be considered a peak.
#' @param min_elution_length Run filter parameter. The minimum length of elution time for a series of signals grouped by 
#'  m/z to be considered a peak.
#' @param mz_tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of 
#'  the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is 
#'  the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline.correct After grouping the observations, the highest intensity in each group is found. If 
#'  the highest is lower than this value, the entire group will be deleted. The default value is NA, in which 
#'  case the program uses the 75th percentile of the height of the noise groups.
#' @param baseline.correct.noise.percentile The perenctile of signal strength of those EIC that don't pass the 
#'  run filter, to be used as the baseline threshold of signal strength.
#' @param do.plot Indicates whether plot should be drawn.
#' @param intensity.weighted Whether to use intensity to weight mass density estimation.
#' @param cache Whether to use cache
#' @return A matrix with four columns: m/z value, retention time, intensity, and group number.
#' @export
#' @examples
#' proc.cdf(input_path, min_pres, min_run, tol, intensity.weighted = intensity_weighted)
proc.cdf <- function(filename,
                     min_presence = 0.5,
                     min_elution_length = 12,
                     mz_tol = 1e-05,
                     baseline_correct = 0.0,
                     baseline_correct_noise_percentile = 0.05,
                     do.plot = FALSE,
                     intensity_weighted = FALSE,
                     cache = FALSE) {
  raw.prof <- load_data(
    filename,
    cache,
    min_elution_length,
    min_presence,
    mz_tol,
    baseline_correct,
    intensity_weighted
  )

  newprof <- cbind(
    raw.prof$features$mz,
    raw.prof$features$rt,
    raw.prof$features$intensities,
    raw.prof$features$grps
  )
  run.sel <- raw.prof$height.rec[which(raw.prof$height.rec[, 2] >= raw.prof$min.count.run * min_presence & raw.prof$height.rec[, 3] > baseline_correct), 1]

  newprof <- newprof[newprof[, 4] %in% run.sel, ]
  new.prof <- cont.index(
    newprof,
    min.pres = min_presence,
    min.run = min_elution_length
  )

  if (do.plot) {
    plot_raw_profile_histogram(
      raw.prof,
      min_presence,
      baseline_correct,
      baseline_correct_noise_percentile,
      mz_tol,
      new.prof
    )
  }

  new_rec_tibble <- tibble::tibble(
    mz = new.prof$new.rec[, 1],
    rt = new.prof$new.rec[, 2],
    intensity = new.prof$new.rec[, 3],
    group_number = new.prof$new.rec[, 4]
  )

  return(new_rec_tibble)
}
