#' Load raw data from file
#' @export
load_file <- function(filename) {
  this <- load.lcms(filename)

  # this could eventually be replaced using drop_na
  na.sel <- c(which(is.na(this$masses)), which(is.na(this$labels)), which(is.na(this$intensi)))
  if (length(na.sel) > 0) {
    na.sel <- unique(na.sel)
    this$masses <- this$masses[-na.sel]
    this$labels <- this$labels[-na.sel]
    this$intensi <- this$intensi[-na.sel]

    warning("there are NA values in the m/z or intensity. Check the file:", filename)
  }
  # TODO
  # this <- tidyr::drop_na(as.data.frame(this))
  return(this)
}

#' Load data either from cache or load raw file and detect peaks.
#' @export
load_data <- function(filename,
                      cache,
                      min.run,
                      min.pres,
                      tol,
                      baseline.correct,
                      intensity.weighted) {
  rawprof_filename <- paste(strsplit(tolower(filename), "\\.")[[1]][1], "_", min.run, "_", min.pres, "_", tol, ".rawprof", sep = "")

  if (cache && file.exists(rawprof_filename)) {
    load(rawprof_filename)
  } else {
    raw.data <- load_file(filename)
    raw.prof <- adaptive.bin(raw.data, min.run = min.run, min.pres = min.pres, tol = tol, baseline.correct = baseline.correct, weighted = intensity.weighted)
  }

  if (cache && !file.exists(rawprof_filename)) {
    save(raw.prof, file = rawprof_filename)
  }

  return(raw.prof)
}

#' Process file and return the profile.
#' @export
#' @examples
#' proc.cdf(input_path, min_pres, min_run, tol, intensity.weighted = intensity_weighted)
proc.cdf <- function(filename,
                     min.pres = 0.5,
                     min.run = 12,
                     tol = 1e-05,
                     baseline.correct = 0.0,
                     baseline.correct.noise.percentile = 0.05,
                     do.plot = FALSE,
                     intensity.weighted = FALSE,
                     cache = FALSE) {
  raw.prof <- load_data(filename, cache, min.run, min.pres, tol, baseline.correct, intensity.weighted)

  newprof <- cbind(raw.prof$masses, raw.prof$labels, raw.prof$intensi, raw.prof$grps)
  run.sel <- raw.prof$height.rec[which(raw.prof$height.rec[, 2] >= raw.prof$min.count.run * min.pres & raw.prof$height.rec[, 3] > baseline.correct), 1]

  newprof <- newprof[newprof[, 4] %in% run.sel, ]
  new.prof <- cont.index(newprof, min.pres = min.pres, min.run = min.run)

  if (do.plot) {
    plot_raw_profile_histogram(
      raw.prof,
      min.pres,
      baseline.correct,
      baseline.correct.noise.percentile,
      tol,
      new.prof
    )
  }

  return(new.prof$new.rec)
}
