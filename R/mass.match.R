#' Compute matches between mz array and specific mass value with a tolerance.
#' @param sample_mz The mz array for which to compute the matching.
#' @param known_mz The mz value with which to match.
#' @param match_tol_ppm Matching tolerance in ppm.
#' @return Indicies of m/z values within the tolerance of any known m/z.
#' @export
#' @examples
#' find_mz_match(
#'  sample_mz = c(10, 20, 21),
#'  known_mz = 20
#' )
find_mz_match <- function(sample_mz, known_mz, match_tol_ppm = 5) {
  matched_mz_idx <- rep(0, nrow(sample_mz))
  match_tol_ppm <- match_tol_ppm / 1e6

  for (i in seq_along(sample_mz)) {
    rel_diff <- abs((sample_mz[i] - known_mz) / sample_mz[i])
    if (min(rel_diff) < match_tol_ppm) {
      matched_mz_idx[i] <- 1
    }
  }
  return(which(matched_mz_idx == 1))
}
