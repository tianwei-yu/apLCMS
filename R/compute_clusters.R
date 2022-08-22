
#' Compute clusters of mz and rt and assign cluster id to individual features.
#' 
#' @description
#' Uses tolerances to group features with mz and rt within the tolerance into clusters,
#' creating larger features from raw data points. Custom tolerances for mz and rt are
#' computed based on the given parameters.
#' @param feature_tables list List of tibbles containing features.
#' Each tibble should have columns [mz, rt, area].
#' @param mz_tol_relative float Relative mz tolerance to use for grouping features.
#' @param mz_tol_absolute float Absolute mz tolerance to use for grouping features.
#' @param mz_max_diff float Maximum difference between featuure mz values to belong to the same cluster.
#' @param rt_tol_relative float Relative retention time tolerance to use for grouping features.
#' @param do.plot bool Plot graphics or not.
#' @return Returns a list with following items:
#' \itemize{
#'   \item feature_tables - list - Feature tables with added columns [sample_id, cluster].
#'   \item rt_tol_relative - float - Newly determined relative rt tolerance.
#'   \item mz_tol_relative - float - Newly determined relative mz tolerance.
#'}
compute_clusters <- function(feature_tables,
                             mz_tol_relative,
                             mz_tol_absolute,
                             mz_max_diff,
                             rt_tol_relative,
                             do.plot = FALSE) {
  number_of_samples <- length(feature_tables)
  all <- concatenate_feature_tables(feature_tables, "rt")

  if (is.na(mz_tol_relative)) {
    mz_tol_relative <- find.tol(
      all$mz,
      mz_max_diff = mz_max_diff,
      do.plot = do.plot
    )
    if (length(mz_tol_relative) == 0) {
      mz_tol_relative <- 1e-5
      warning("Automatic tolerance finding failed, 10 ppm was assigned.
                        May need to manually assign alignment mz tolerance level.")
    }
  } else if (do.plot) {
    draw_plot(
      main = "m/z tolerance level given",
      label = mz_tol_relative
    )
  }

  if (!is.na(rt_tol_relative) && do.plot) {
    draw_plot(
      main = "retention time \n tolerance level given",
      label = rt_tol_relative
    )
  }

  res <- find.tol.time(
    all,
    number_of_samples = number_of_samples,
    mz_tol_relative = mz_tol_relative,
    rt_tol_relative = rt_tol_relative,
    mz_tol_absolute = mz_tol_absolute,
    do.plot = do.plot
  )
  all.ft <- res$features
  rt_tol_relative <- res$rt.tol

  message("**** performing time correction ****")
  message(paste("m/z tolerance level: ", mz_tol_relative))
  message(paste("time tolerance level:", rt_tol_relative))

  # Select features from individual samples, sort by mz and rt and 
  # return the sorted tables as individual tibbles.
  feature_tables <- res$features |>
    dplyr::group_by(sample_id) |>
    dplyr::arrange_at(c("mz", "rt")) |>
    dplyr::group_split()

  return(list(feature_tables = feature_tables, rt_tol_relative = rt_tol_relative, mz_tol_relative = mz_tol_relative))
}