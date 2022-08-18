#' @import dplyr stringr tibble tools tidyr parallel doParallel snow
NULL
#> NULL

merge_known_tables <- function(batchwise, batches_idx) {
  colnames <- c("chemical_formula", "HMDB_ID", "KEGG_compound_ID", "mass", "ion.type", "m.z",
              "Number_profiles_processed", "Percent_found", "mz_min", "mz_max", 
              "RT_mean", "RT_sd", "RT_min", "RT_max", "int_mean(log)", "int_sd(log)", 
              "int_min(log)", "int_max(log)")

  known_table <- tibble(
    chemical_formula = character(),
    HMDB_ID = character(),
    KEGG_compound_ID = character(),
    mass = numeric(),
    ion.type = character(),
    m.z = numeric(),
    Number_profiles_processed = numeric(),
    Percent_found = numeric(),
    mz_min = numeric(),
    mz_max = numeric(),
    RT_mean = numeric(),
    RT_sd = numeric(),
    RT_min = numeric(),
    RT_max = numeric(),
    "int_mean(log)" = numeric(),
    "int_sd(log)" = numeric(),
    "int_min(log)" = numeric(),
    "int_max(log)" = numeric()
  )

  for (batch in batches_idx) {
    known_table <- dplyr::full_join(known_table, batchwise[[batch]]$updated.known.table, by = colnames)
  }

  return(known_table)
}

filter_features_by_presence <- function(feature_table,
                                        metadata,
                                        batches_idx,
                                        within_batch_threshold,
                                        across_batch_threshold) {
  across_batch_presence <- data.frame(
    matrix(
      nrow = nrow(feature_table),
      ncol = length(batches_idx)
    )
  )
  for (batch_id in batches_idx) {
    samples <- filter(metadata, batch == batch_id)$sample_name
    intensities <- dplyr::select(feature_table, contains("_intensity"))
    presence <- intensities > 0
    above_threshold <- rowSums(presence) / length(samples) >= within_batch_threshold
    across_batch_presence[, batch_id] <- above_threshold
  }
  above_threshold <- rowSums(across_batch_presence) / length(batches_idx) >= across_batch_threshold
  return(feature_table[above_threshold, ])
}

readjust_times <- function(within_batch, between_batch) {
  within_batch_recovered <- long_to_wide_feature_table(
    within_batch$recovered_feature_sample_table
  )
  between_batch_rts <- between_batch$rt
  for (j in 1:length(within_batch$corrected_features)) {
    for (i in 1:nrow(within_batch$corrected_features[[j]])) {
      diff.time <- abs(
        within_batch_recovered$rt -
          within_batch$corrected_features[[j]][i, "pos"]
      )
      min_idx <- which(diff.time == min(diff.time))[1]
      within_batch$corrected_features[[j]][i, "pos"] <- between_batch_rts[min_idx]
    }
  }
  return(within_batch$corrected_features)
}

compute_intensity_medians <- function(feature_table) {
  stopifnot("sample_intensity" %in% colnames(feature_table))
  feature_table <- dplyr::group_by(feature_table, mz, rt) %>%
    dplyr::mutate(median_intensity = median(sample_intensity)) %>%
    dplyr::ungroup()
  return(feature_table)
}

bind_batch_label_column <- function(filenames, metadata) {
  stopifnot(nrow(metadata) == length(filenames))

  filenames <- as_tibble(filenames)
  colnames(filenames)[1] <- "filename"
  filenames <- mutate(filenames, sample_name = get_sample_name(filename))
  filenames <- inner_join(filenames, metadata, by = "sample_name")
  return(dplyr::select(filenames, -sample_name))
}

feature_recovery <- function(cluster,
                             step_one_features,
                             batchwise,
                             filenames_batchwise,
                             corrected,
                             aligned,
                             batches_idx,
                             mz.tol,
                             batch.align.mz.tol,
                             batch.align.rt.tol,
                             recover.mz.range,
                             recover.rt.range,
                             use.observed.range,
                             min.bw,
                             max.bw,
                             recover.min.count) {
  recovered <- tibble(
    mz = numeric(),
    rt = numeric(),
    mz_min = numeric(),
    mz_max = numeric()
  )

  for (batch_id in batches_idx)
  {
    this.fake <- long_to_wide_feature_table(step_one_features[[batch_id]])
    this.fake.medians <- apply(dplyr::select(this.fake, contains("_intensity")), 1, median)

    # adjusting the time (already within batch adjusted)
    this.features <- readjust_times(batchwise[[batch_id]], corrected[[batch_id]])
    aligned_intensities <- dplyr::select(aligned, contains("_intensity"))
    batchwise_intensities <- dplyr::select(this.fake, contains("_intensity"))

    this.fake.time <- dplyr::select(this.fake, contains("_rt"))
    this.pk.time <- this.aligned <- matrix(0, nrow = nrow(aligned), ncol = ncol(batchwise_intensities))

    for (sample in 1:nrow(aligned)) {
      if (aligned_intensities[sample, batch_id] != 0) {
        idx <- which(between(this.fake$mz, aligned[sample, "mz_min"], aligned[sample, "mz_max"]) &
          abs(this.fake.medians - aligned_intensities[sample, batch_id]) < 1)
        if (length(idx) < 1) {
          idx <- which(between(this.fake$mz, aligned[sample, "mz_min"], aligned[sample, "mz_max"]))
        }
        if (length(idx) < 1) {
          message("warning: batch ", batch_id, " sample ", sample, " has matching issue")
        } else {
          this.aligned[sample, ] <- apply(batchwise_intensities[idx, ], 2, sum)
          this.pk.time[sample, ] <- apply(this.fake.time[idx, ], 2, median)
        }
      } else {
        ### go into individual feature tables to find a match
        recaptured <- rep(0, ncol(this.aligned))
        recaptured.time <- rep(NA, ncol(this.aligned))

        for (j in 1:length(this.features)) {
          diff.mz <- abs(this.features[[j]][, "mz"] - aligned[sample, "mz"])
          diff.time <- abs(this.features[[j]][, "pos"] - aligned[sample, "rt"])
          idx <- which(diff.mz < aligned[sample, "mz"] * batch.align.mz.tol & diff.time <= batch.align.rt.tol)

          if (length(idx) > 0) {
            idx <- idx[which(diff.time[idx] == min(diff.time[idx]))[1]]
            recaptured[j] <- this.features[[j]][idx, "area"]
            recaptured.time[j] <- this.features[[j]][idx, "pos"]
          }
        }
        this.aligned[sample, ] <- recaptured
        this.pk.time[sample, ] <- recaptured.time
      }
    }

    colnames(this.aligned) <- extract_pattern_colnames(this.fake, "_intensity")
    colnames(this.aligned) <- stringr::str_remove_all(colnames(this.aligned), "_intensity")
    colnames(this.pk.time) <- extract_pattern_colnames(this.fake, "_rt")
    colnames(this.pk.time) <- stringr::str_remove_all(colnames(this.pk.time), "_rt")


    aligned_features <- dplyr::select(aligned, mz, rt, mz_min, mz_max)
    aligned_int_crosstab <- dplyr::bind_cols(aligned_features, as_tibble(this.aligned))
    aligned_rt_crosstab <- dplyr::bind_cols(aligned_features, as_tibble(this.pk.time))
    recovered_batchwise <- recover_weaker_signals(
      cluster = cluster,
      filenames = filter(filenames_batchwise, batch == batch_id)$filename,
      extracted_features = batchwise[[batch_id]]$extracted_features,
      corrected_features = batchwise[[batch_id]]$corrected_features,
      aligned_rt_crosstab = aligned_rt_crosstab,
      aligned_int_crosstab = aligned_int_crosstab,
      original_mz_tolerance = mz.tol,
      aligned_mz_tolerance = batch.align.mz.tol,
      aligned_rt_tolerance = batch.align.rt.tol,
      mz_range = recover.mz.range,
      rt_range = recover.rt.range,
      use_observed_range = use.observed.range,
      min_bandwidth = min.bw,
      max_bandwidth = max.bw,
      recover_min_count = recover.min.count
    )

    recovered_batchwise <- as_wide_aligned_table(recovered_batchwise)

    recovered <- dplyr::full_join(recovered, recovered_batchwise, by = c("mz", "rt", "mz_min", "mz_max"))
  }
  recovered <- dplyr::select(recovered, mz, rt, mz_min, mz_max, contains("_rt"), contains("_intensity"))

  return(recovered)
}

semisup_to_hybrid_adapter <- function(batchwise, batches_idx) {
  for (batch in batches_idx) {
    final.ftrs <- as_tibble(batchwise[[batch]]$final.ftrs)
    final.times <- as_tibble(batchwise[[batch]]$final.times)

    mz_pattern <- c("mz.min", "mz.max")
    mz_replacement <- c("mz_min", "mz_max")

    colnames(final.ftrs) <- stringr::str_replace_all(
      colnames(final.ftrs),
      mz_pattern,
      mz_replacement)

    colnames(final.times) <- stringr::str_replace_all(
      colnames(final.ftrs),
      mz_pattern,
      mz_replacement )

    feature_cols <- c("mz", "time", "mz_min", "mz_max")
    sample_cols_idx <- which(!colnames(final.ftrs) %in% feature_cols)

    colnames(final.ftrs)[sample_cols_idx] <- tools::file_path_sans_ext(colnames(final.ftrs)[sample_cols_idx])
    colnames(final.times)[sample_cols_idx] <- tools::file_path_sans_ext(colnames(final.times)[sample_cols_idx])

    sample_cols <- colnames(final.ftrs)[sample_cols_idx]

    final.ftrs <- tidyr::pivot_longer(final.ftrs,
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "sample_intensity"
    )

    final.times <- tidyr::pivot_longer(final.times,
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "sample_rt"
    )

    recovered_feature_sample_table <- dplyr::full_join(final.times, final.ftrs,
      by = c("mz", "time", "mz_min", "mz_max", "sample")
    )

    colnames(recovered_feature_sample_table)[2] <- "rt"

    batchwise[[batch]]$recovered_feature_sample_table <- recovered_feature_sample_table

    batchwise[[batch]]$corrected_features <- batchwise[[batch]]$features2
    batchwise[[batch]]$extracted_features <- batchwise[[batch]]$features
  }

  return(batchwise)
}

#' Two step hybrid feature detection.
#' 
#' A two-stage hybrid feature detection and alignment procedure, for data generated in multiple batches.
#' 
#' @param filenames file names
#' @param metadata the batch label of each file.
#' @param work_dir The folder where all CDF files to be processed are located.
#' @param min.within.batch.prop.detect A feature needs to be present in at least this proportion of the files, 
#'  for it to be initially detected as a feature for a batch. This parameter replaces the "min.exp" parameter in semi.sup().
#' @param min.within.batch.prop.report A feature needs to be present in at least this proportion of the files, 
#'  in a proportion of batches controlled by "min.batch.prop", to be included in the final feature table. This parameter 
#'  replaces the "min.exp" parameter in semi.sup().
#' @param min.batch.prop A feature needs to be present in at least this proportion of the batches, for it to be 
#'  considered in the entire data.
#' @param batch.align.mz.tol The m/z tolerance in ppm for between-batch alignment.
#' @param batch.align.rt.tol The RT tolerance for between-batch alignment.
#' @param known.table A data frame containing the known metabolite ions and previously found features.
#' @param cluster The number of CPU cores to be used
#' @param min.pres This is a parameter of the run filter, to be passed to the function proc.cdf().
#' @param min.run This is a parameter of the run filter, to be passed to the function proc.cdf().
#' @param mz.tol The user can provide the m/z tolerance level for peak identification. This value is expressed as the 
#'  percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param baseline.correct.noise.percentile The perenctile of signal strength of those EIC that don't pass the run filter, 
#'  to be used as the baseline threshold of signal strength. This parameter is passed to proc.cdf()
#' @param shape.model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian". 
#'  When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param baseline.correct This is a parameter in peak detection. After grouping the observations, the highest observation 
#'  in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, 
#'  which allows the program to search for the cutoff level.
#' @param peak.estim.method the bi-Gaussian peak parameter estimation method, to be passed to subroutine prof.to.features. 
#'  Two possible values: moment and EM.
#' @param min.bw The minimum bandwidth in the smoother in prof.to.features().
#' @param max.bw The maximum bandwidth in the smoother in prof.to.features().
#' @param sd.cut A parameter for the prof.to.features() function. A vector of two. Features with standard deviation outside 
#'  the range defined by the two numbers are eliminated.
#' @param sigma.ratio.lim A parameter for the prof.to.features() function. A vector of two. It enforces the belief of the 
#'  range of the ratio between the left-standard deviation and the right-standard deviation of the bi-Gaussian function used 
#'  to fit the data.
#' @param component.eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts 
#'  for a proportion of intensities less than this value, the component will be ignored.
#' @param moment.power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param align.mz.tol The user can provide the m/z tolerance level for peak alignment to override the program's selection. 
#'  This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param align.rt.tol The user can provide the elution time tolerance level to override the program's selection. This value 
#'  is in the same unit as the elution time, normaly seconds.
#' @param max.align.mz.diff As the m/z tolerance in alignment is expressed in relative terms (ppm), it may not be suitable 
#'  when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly influences feature matching 
#'  in higher m/z range.
#' @param pre.process Logical. If true, the program will not perform time correction and alignment. It will only generate peak tables 
#'  for each spectra and save the files. It allows manually dividing the task to multiple machines.
#' @param recover.mz.range A parameter of the recover.weaker() function. The m/z around the feature m/z to search for observations. 
#'  The default value is NA, in which case 1.5 times the m/z tolerance in the aligned object will be used.
#' @param recover.rt.range A parameter of the recover.weaker() function. The retention time around the feature retention time to 
#'  search for observations. The default value is NA, in which case 0.5 times the retention time tolerance in the aligned 
#'  object will be used.
#' @param use.observed.range A parameter of the recover.weaker() function. If the value is TRUE, the actual range of the observed 
#'  locations of the feature in all the spectra will be used.
#' @param match.tol.ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param new.feature.min.count The number of profiles a new feature must be present for it to be added to the database.
#' @param recover.min.count The minimum time point count for a series of point in the EIC for it to be considered a true feature.
#' @param intensity.weighted Whether to use intensity to weight mass density estimation.
#' @param BIC.factor the factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1, 
#'  models with more peaks are penalized more.
#' @return A list is returned.
#' \itemize{
#'   \item batchwise.results - A list. Each item in the list is the product of semi.sup() from a single batch.
#'   \item final.ftrs - Feature table. This is the end product of the function.
#' }
#' @export
#' @examples
#' two.step.hybrid(test_names, metadata, tempdir, known.table = known_table, cluster = num_workers)
two.step.hybrid <- function(filenames,
                            metadata,
                            work_dir,
                            min.within.batch.prop.detect = 0.1,
                            min.within.batch.prop.report = 0.5,
                            min.batch.prop = 0.5,
                            batch.align.mz.tol = 1e-5,
                            batch.align.rt.tol = 50,
                            known.table = NA,
                            cluster = 4,
                            min.pres = 0.5,
                            min.run = 12,
                            mz.tol = 1e-5,
                            baseline.correct.noise.percentile = 0.05,
                            shape.model = "bi-Gaussian",
                            baseline.correct = 0,
                            peak.estim.method = "moment",
                            min.bw = NA,
                            max.bw = NA,
                            sd.cut = c(0.1, 100),
                            sigma.ratio.lim = c(0.05, 20),
                            component.eliminate = 0.01,
                            moment.power = 2,
                            align.mz.tol = NA,
                            align.rt.tol = NA,
                            max.align.mz.diff = 0.01,
                            pre.process = FALSE,
                            recover.mz.range = NA,
                            recover.rt.range = NA,
                            use.observed.range = TRUE,
                            match.tol.ppm = NA,
                            new.feature.min.count = 2,
                            recover.min.count = 3,
                            intensity.weighted = FALSE,
                            BIC.factor = 2) {
  filenames_batchwise <- bind_batch_label_column(filenames, metadata)
  batches_idx <- unique(metadata$batch)
  batchwise <- new("list")
  message("* processing ", length(batches_idx), " batches separately")

  for (batch.i in batches_idx) {
    files_batch <- dplyr::filter(filenames_batchwise, batch == batch.i)$filename
    message("* processing ", length(files_batch), " samples from batch ", batch.i)

    features <- semi.sup(
      files = files_batch,
      folder = work_dir,
      n.nodes = cluster,
      known.table = known.table,
      sd.cut = sd.cut,
      sigma.ratio.lim = sigma.ratio.lim,
      component.eliminate = component.eliminate,
      moment.power = moment.power,
      min.pres = min.pres,
      min.run = min.run,
      min.exp = ceiling(min.within.batch.prop.detect * length(files_batch)),
      mz.tol = mz.tol,
      baseline.correct.noise.percentile = baseline.correct.noise.percentile,
      baseline.correct = baseline.correct,
      align.mz.tol = align.mz.tol,
      align.rt.tol = align.rt.tol,
      max.align.mz.diff = max.align.mz.diff,
      recover.mz.range = recover.mz.range,
      recover.rt.range = recover.rt.range,
      use.observed.range = use.observed.range,
      shape.model = shape.model,
      new.feature.min.count = new.feature.min.count,
      recover.min.count = recover.min.count
    )

    features$final.ftrs <- features$final.ftrs[order(features$final.ftrs[, 1], features$final.ftrs[, 2]), ]
    features$final.times <- features$final.times[order(features$final.times[, 1], features$final.times[, 2]), ]

    batchwise[[batch.i]] <- features
  }

  batchwise <- semisup_to_hybrid_adapter(batchwise, batches_idx)

  step_one_features <- list()
  for (batch_id in batches_idx) {
    step_one_features[[batch_id]] <- compute_intensity_medians(
      batchwise[[batch_id]]$recovered_feature_sample_table
    ) %>%
      dplyr::select(-c("mz_min", "mz_max"))
  }

  cluster <- parallel::makeCluster(cluster)
  doParallel::registerDoParallel(cluster)
  snow::clusterEvalQ(cluster, library(recetox.aplcms))

  extracted_features <- list()
  for (batch_id in batches_idx) {
    extracted <- batchwise[[batch_id]]$recovered_feature_sample_table
    extracted <- long_to_wide_feature_table(extracted)
    intensities <- dplyr::select(extracted, contains("_intensity"))
    median_intensity <- apply(intensities, 1, median)
    extracted <- dplyr::select(extracted, mz, rt)
    extracted[, 3:4] <- NA
    colnames(extracted)[3:4] <- c("mz_min", "mz_max")
    extracted[, 5] <- median_intensity
    colnames(extracted)[5] <- "median_intensity"
    extracted_features[[batch_id]] <- extracted
  }

  message("* aligning time")
  corrected <- adjust.time(extracted_features,
    mz_tol_relative = batch.align.mz.tol,
    rt_tol_relative = batch.align.rt.tol,
    mz_max_diff = 10 * mz.tol,
    mz_tol_absolute = max.align.mz.diff,
    rt_colname = "rt"
  )

  message("* aligning features")
  aligned <- align_features(
    sample_names = paste0("batch_", batches_idx),
    features = corrected,
    min_occurrence = ceiling(min.batch.prop * length(batches_idx)),
    mz_tol_relative = batch.align.mz.tol,
    rt_tol_relative = batch.align.rt.tol,
    mz_max_diff = 10 * mz_tol,
    mz_tol_absolute = max.align.mz.diff,
    rt_colname = "rt"
  )

  aligned_wide <- as_wide_aligned_table(aligned)

  message("* recovering features across batches")
  recovered <- feature_recovery(
    cluster = cluster,
    step_one_features = step_one_features,
    batchwise = batchwise,
    filenames_batchwise = filenames_batchwise,
    corrected = corrected,
    aligned = aligned_wide,
    batches_idx = batches_idx,
    mz.tol = mz.tol,
    batch.align.mz.tol = batch.align.mz.tol,
    batch.align.rt.tol = batch.align.rt.tol,
    recover.mz.range = recover.mz.range,
    recover.rt.range = recover.rt.range,
    use.observed.range = use.observed.range,
    min.bw = min.bw,
    max.bw = max.bw,
    recover.min.count = recover.min.count
  )

  snow::stopCluster(cluster)

  recovered_filtered <- filter_features_by_presence(
    feature_table = recovered,
    metadata = metadata,
    batches_idx = batches_idx,
    within_batch_threshold = min.within.batch.prop.report,
    across_batch_threshold = min.batch.prop
  )

  recovered_filtered <- dplyr::arrange(recovered_filtered, mz, rt)
  recovered_features <- wide_to_long_feature_table(recovered_filtered)

  features <- new("list")
  features$batchwise_features <- batchwise
  features$known_table <- merge_known_tables(batchwise, batches_idx)
  features$aligned_features <- as_feature_sample_table(
    rt_crosstab = aligned$rt_crosstab,
    int_crosstab = aligned$int_crosstab
  )
  features$corrected_features <- corrected
  features$final_features <- recovered_features
  return(features)
}
