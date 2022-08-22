#' @import mzR
#' @import tibble
NULL
#> NULL

#' Loading LC/MS data.
#' 
#' This is an internal function. It loads LC/MS data into memory.
#' 
#' @param filename The CDF file name.
#' @return A list is returned.
#' \itemize{
#'   \item masses - The vector of m/z values.
#'   \item labels - The vector of retention times.
#'   \item intensi - The vector of intensity values.
#'   \item times - The vector of unique time points.
#' }
#' @export
#' @examples
#' load.lcms(filename)
load.lcms <- function(filename) {
  mz_conn <- mzR::openMSfile(filename = filename)
  b <- mzR::header(mz_conn)$retentionTime

  masses <- NULL
  intensi <- NULL
  labels <- NULL

  segs <- seq(0, length(b), by = 200)
  if ((length(b) %% 200) != 0)
    segs <- c(segs, length(b))

  for (n in 2:length(segs)) {
    a <- mzR::peaks(mz_conn, scans = (segs[n - 1] + 1):segs[n])

    this_masses <- NULL
    this_intensi <- NULL
    this_labels <- NULL

    for (i in seq_along(a)) {
      this_a <- a[[i]]

      if (!is.null(nrow(this_a))) {
        this_a <- this_a[this_a[, 2] > 1e-10,]
        if (is.null(nrow(this_a)))
          this_a <- matrix(this_a, nrow = 1)

        this_masses <- c(this_masses, this_a[, 1])
        this_intensi <- c(this_intensi, this_a[, 2])
        this_labels <- c(this_labels, rep(b[segs[n - 1] + i], nrow(this_a)))
      }else {
        b[segs[n - 1] + i] <- NA
      }
    }

    masses <- c(masses, this_masses)
    intensi <- c(intensi, this_intensi)
    labels <- c(labels, this_labels)
  }

  mzR::close(mz_conn)

  features <- tibble::tibble(mz = masses, rt = labels, intensities = intensi)

  return(features)
}
