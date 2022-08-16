#' @import dplyr tidyr tibble stringr arrow
NULL
#> NULL

#' Concatenate multiple feature lists and add the sample id (origin of feature) as additional column.
#' 
#' @param features list List of tibbles containing extracted feature tables.
#' @param rt_colname string Name of retention time information column, usually "pos".
concatenate_feature_tables <- function(features, rt_colname) {
  for (i in seq_along(features)) {
      if(!("sample_id" %in% colnames(features[[i]]))) {
        features[[i]] <- tibble::add_column(features[[i]], sample_id = i)
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

load_aligned_features <- function(rt_file, int_file, tol_file) {
  rt_cross_table <- arrow::read_parquet(rt_file)
  int_cross_table <- arrow::read_parquet(int_file)
  tolerances_table <- arrow::read_parquet(tol_file)
  
  row.names(rt_cross_table) <- as.character(row.names(rt_cross_table))
  row.names(int_cross_table) <- as.character(row.names(int_cross_table))
  
  result <- list()
  result$mz_tolerance <- tolerances_table$mz_tolerance
  result$rt_tolerance <- tolerances_table$rt_tolerance
  result$rt_crosstab <- rt_cross_table
  result$int_crosstab <- int_cross_table
  return(result)
}

create_feature_sample_table <- function(features) {
  table <- as_feature_sample_table(
    rt_crosstab = features$rt_crosstab,
    int_crosstab = features$int_crosstab
  )
  return(table)
}

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