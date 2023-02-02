#' @import dplyr tidyr tibble stringr arrow
NULL
#> NULL

register_functions_to_cluster <- function(cluster) {
    snow::clusterExport(cluster, list(
        'proc.cdf',
        'prof.to.features',
        'load.lcms',
        'adaptive.bin',
        'find.turn.point',
        'msExtrema',
        'find_local_maxima',
        'combine.seq.3',
        'run_filter',
        'interpol.area',
        'load_file',
        'load_data',
        'plot_raw_profile_histogram',
        'compute_mass_values',
        'compute_densities',
        'compute_breaks',
        'compute_breaks_3',
        'compute_boundaries',
        'compute_rt_intervals_indices',
        'increment_counter',
        'rm.ridge',
        'compute_delta_rt',
        'bigauss.mix',
        'bigauss.esti',
        'rev_cum_sum',
        'compute_bounds',
        'validate_model_method_input',
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
        'compute_scale',
        'span',
        'compute_uniq_grp',
        'predict_smoothed_rt',
        'label_val_to_keep'
    ))
    snow::clusterEvalQ(cluster, library("dplyr"))
}

#' Concatenate multiple feature lists and add the sample id (origin of feature) as additional column.
#' 
#' @param features list List of tibbles containing extracted feature tables.
concatenate_feature_tables <- function(features, sample_names) {
    for (i in seq_along(features)) {
        if(!("sample_id" %in% colnames(features[[i]]))) {
            features[[i]] <- tibble::add_column(features[[i]], sample_id = sample_names[i])
        }
    }
    
    merged <- dplyr::bind_rows(features)
    return(merged)
}

#' @export
load_aligned_features <- function(metadata_file, intensities_file, rt_file, tol_file) {
    metadata <- arrow::read_parquet(metadata_file)
    intensities <- arrow::read_parquet(intensities_file)
    rt <- arrow::read_parquet(rt_file)
    tolerances <- arrow::read_parquet(tol_file)
    
    result <- list()
    result$metadata <- as_tibble(metadata)
    result$intensity <- as_tibble(intensities)
    result$rt <- as_tibble(rt)
    result$mz_tol_relative <- tolerances$mz_tolerance
    result$rt_tol_relative <- tolerances$rt_tolerance
    return(result)
}

#' @export
span <- function(x) {
    diff(range(x, na.rm = TRUE))
}

#' @description
#' Compute standard deviation of m/z values groupwise
#' @export
compute_mz_sd <- function(feature_groups) {
    mz_sd <- c()
    for (i in seq_along(feature_groups)) {
        group <- feature_groups[[i]]
        
        if (nrow(group > 1)) {
            group_sd <- sd(group[, "mz"])
            mz_sd <- append(mz_sd, group_sd)
        }
    }
    return(mz_sd)
}

#' @export
get_num_workers <- function() {
    # CRAN limits the number of cores available to packages to 2
    # source https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions#50571533
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    
    if (nzchar(chk) && chk == "TRUE") {
        # use 2 cores in CRAN/Travis/AppVeyor
        num_workers <- 2L
    } else {
        # use all cores in devtools::test()
        num_workers <- parallel::detectCores()
    }
    return(num_workers)
}