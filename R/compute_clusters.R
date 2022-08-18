
add_sample_id_and_rt_cluster <- function(sample, all.ft, current_sample_id) {
  sample <- sample |> dplyr::arrange_at(c("mz", "rt"))
  sample_grouped <- all.ft |> dplyr::filter(sample_id == current_sample_id) |> dplyr::arrange_at(c("mz", "rt"))
  
  if(!tibble::has_name(sample, "sample_id")) {
    sample <- tibble::add_column(sample, sample_id = current_sample_id)
  }

  if(tibble::has_name(sample, "cluster")) {
    sample <- sample |> dplyr::select(-cluster)
  }

  features <- dplyr::bind_cols(sample, dplyr::select(sample_grouped, cluster))
  return(features)
}

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

  for (i in 1:number_of_samples) {
    features <- add_sample_id_and_rt_cluster(
      feature_tables[[i]],
      all.ft,
      i
    )
    feature_tables[[i]] <- features
  }
  return(list(feature_tables = feature_tables, rt_tol_relative = rt_tol_relative, mz_tol_relative = mz_tol_relative))
}