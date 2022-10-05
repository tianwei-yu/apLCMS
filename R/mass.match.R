#' Compute matches between mz array and specific mass value with a tolerance.
#' @param sample_mz The mz array for which to compute the matching.
#' @param known_mz The mz value with which to match.
#' @param match_tol_ppm Matching tolerance in ppm.
#' @return Binary vector, 1 indicating a match, 0 a mismatch.
#' @export
#' @examples
#' mass.match(
#'  sample_mz = c(10, 20, 21),
#'  known_mz = 20
#' )
mass.match <- function(sample_mz, known_mz, match_tol_ppm = 5) {
  matched_mz_idx <- rep(0, nrow(sample_mz))
  for (i in seq_along(sample_mz))
  {
    this.d <- abs((sample_mz[i] - known_mz) / sample_mz[i])
    if (min(this.d) < match_tol_ppm / 1e6) {
      matched_mz_idx[i] <- 1
    }
  }
  return(matched_mz_idx)
}
