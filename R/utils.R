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
        'compute_scale',
        'span'
    ))
    snow::clusterEvalQ(cluster, library("dplyr"))
}


#' @import dplyr tidyr tibble stringr arrow
NULL
#> NULL

#' Concatenate multiple feature lists and add the sample id (origin of feature) as additional column.
#' 
#' @param features list List of tibbles containing extracted feature tables.
concatenate_feature_tables <- function(features, sample_names) {
    if(!all(is.na(sample_names))) {
        for (i in seq_along(features)) {
            if(!("sample_id" %in% colnames(features[[i]]))) {
                features[[i]] <- tibble::add_column(features[[i]], sample_id = sample_names[i])
            }
        }
    }
    
    merged <- dplyr::bind_rows(features)
    return(merged)
}

#' @export
extract_pattern_colnames <- function(dataframe, pattern) {
    dataframe <- dplyr::select(dataframe, contains(pattern))
    return(colnames(dataframe))
}

as_wide_aligned_table <- function(aligned) {
    mz_scale_table <- aligned$rt_crosstab[, c("mz", "rt", "mz_min", "mz_max")]
    aligned <- as_feature_sample_table(
        rt_crosstab = aligned$rt_crosstab,
        int_crosstab = aligned$int_crosstab
    )
    aligned <- long_to_wide_feature_table(aligned)
    aligned <- dplyr::inner_join(aligned, mz_scale_table, by = c("mz", "rt"))
    return(aligned)
}

pivot_feature_values <- function(feature_table, variable) {
    extended_variable <- paste0("sample_", variable)
    values <- dplyr::select(feature_table, mz, rt, sample, !!sym(extended_variable))
    values <- tidyr::pivot_wider(values, names_from = sample, values_from = !!sym(extended_variable))
    variable_colnames <- colnames(values)[3:ncol(values)]
    variable_colnames <- paste0(variable_colnames, "_", variable)
    colnames(values)[3:ncol(values)] <- variable_colnames
    return(values)
}

long_to_wide_feature_table <- function(feature_table) {
    sample_rts <- pivot_feature_values(feature_table, "rt")
    sample_intensities <- pivot_feature_values(feature_table, "intensity")
    feature_table <- dplyr::select(feature_table, mz, rt) %>%
        dplyr::distinct(mz, rt) %>%
        dplyr::inner_join(sample_rts, by = c("mz", "rt")) %>%
        dplyr::inner_join(sample_intensities, by = c("mz", "rt"))
}

wide_to_long_feature_table <- function(wide_table, sample_names) {
    wide_table <- tibble::rowid_to_column(wide_table, "feature")
    
    long_rt <- tidyr::gather(wide_table, sample, sample_rt, contains("_rt"), factor_key=FALSE) %>%
        dplyr::select(-contains("_intensity")) %>%
        mutate(sample = stringr::str_remove_all(sample, "_rt"))
    long_int <- tidyr::gather(wide_table, sample, sample_intensity, contains("_intensity"), factor_key=FALSE) %>%
        dplyr::select(-contains("_rt")) %>%
        mutate(sample = stringr::str_remove_all(sample, "_intensity"))
    
    long_features <- dplyr::full_join(long_rt, long_int, by = c("feature", "mz", "rt", "mz_min", "mz_max", "sample"))
    
    return(long_features)
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

create_feature_sample_table <- function(features) {
    table <- as_feature_sample_table(
        rt_crosstab = features$rt,
        int_crosstab = features$intensity
    )
    return(table)
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