#' @description
#' Computes unique groups
#' @param min_count_run filter parameter. 
#' @param min_pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped
#'  by m/z to be considered a peak.
#' @param profile The matrix containing m/z, retention time, intensity, and EIC label as columns.
#' @return unique_grp. 
compute_uniq_grp <- function(profile, min_count_run, min_pres = 0.6) {
  grps <- profile  
  ttt <- table(grps)
  ttt <- ttt[ttt >= max(min_count_run * min_pres, 2)]
  unique_grp <- as.numeric(names(ttt))
  return(unique_grp)
}

#' @description
#' Computes the smoothed retention times by using The Nadaraya-Watson kernel regression estimate function.
#' @param min_run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param times. Retention times vector.
#' @return predicted rt.
#' @examples
#' predict_smoothed_rt(min_run = min_run, times)
predict_smoothed_rt <- function(min_run = 5, times) {
  # ksmooth(x, y, kernel, bandwidth, range, n.points, x.points)
  smooth <- ksmooth(
    seq(-min_run + 1, length(times) + min_run), 
    c(rep(0, min_run),   
    times,       
    rep(0, min_run)),    
    kernel = "box",      
    bandwidth = min_run, 
    x.points = 1:length(times) 
  ) 
  # vector of smoothed estimates for the regression at the corresponding x
  smooth <- smooth$y  
  return(smooth)
}

#' @description
#' This function labels the indices of values kept to perform further calculations
#' @param min_run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param min_pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped
#'  by m/z to be considered a peak.
#' @param timeline. 
#' @param this_times. 
#' @param times. Retention times vector.
#' @return to_keep. 
label_val_to_keep <- function(min_run, timeline, min_pres, this_times, times) {
    this_timeline <- timeline
    this_timeline[this_times] <- 1
    to_keep <- this_times * 0
    
    # filtering based on the kernel regression estimate
    this_smooth <- predict_smoothed_rt(min_run, this_timeline)
    if (max(this_smooth) >= min_pres) {
      measured_points <- good_points <- timeline
      measured_points[this_times] <- 1

      good_sel <- which(this_smooth >= min_pres)
      good_points[good_sel] <- 1
      for (j in (-min_run):min_run) {
        curr_sel <- good_sel + j
        curr_sel <- curr_sel[curr_sel > 0 & curr_sel <= length(times)]
        good_points[curr_sel] <- 1
      }

      measured_points <- measured_points * good_points
      to_keep[which(this_times %in% which(measured_points == 1))] <- 1     
    }
    return(to_keep)
}
#' @description
#' Continuity index. 
#' Internal function that removes noise in the retention time dimension. It uses continuity index (or "run filter") to select putative peaks from EIC. 
#' @param newprof The matrix containing m/z, retention time, intensity, and EIC label as columns.
#' @param min_pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped
#' by m/z to be considered a peak.
#' @param min_run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @return A list is returned. new_rec - The matrix containing m/z, retention time, intensity, and EIC label as columns after applying the run filter.
#' @export
#' @examples
#' run_filter(newprof, min_pres = min_pres, min_run = min_run)
run_filter <- function(newprof,
                       min_pres = 0.6,
                       min_run = 5) {
  
  newprof <- tibble::tibble(mz = newprof[,1], rt = newprof[,2], intensi = newprof[,3], grps = newprof[,4])

  # ordering retention time values
  labels <- newprof$rt
  times <- unique(labels)
  times <- times[order(times)]  
 
  for (i in 1:length(times)) labels[which(newprof$rt == times[i])] <- i #now labels is the index of time points
  newprof$rt <- labels  

  # calculates the minimun number of rt points to be considered a peak
  min_count_run <- min_run * length(times) / (max(times) - min(times))
  min_run <- round(min_count_run)

  # computes unique groups 
  uniq_grp <- compute_uniq_grp(newprof$grps, min_count_run, min_pres)
  
  # ordered by mz and grps data that are inside unigrps
  newprof <- dplyr::filter(newprof, grps %in% uniq_grp) |> dplyr::arrange(grps, mz)
  
  # computes break points i.e. indices of mass differences greater than min_mz_tol
  breaks <- compute_breaks_3(newprof$grps)

  # init counters for loop
  new_rec <- newprof * 0
  rec_pointer <- 1
  timeline <- rep(0, length(times))
  for (m in 2:length(breaks))
  {
    this_prof <- dplyr::slice(newprof, (breaks[m - 1] + 1):breaks[m]) |> dplyr::arrange_at("rt")

    to_keep <- label_val_to_keep(
      min_run,
      timeline,
      min_pres,
      this_prof$rt,
      times
    )
    
    # operation over selected indices 
    if (sum(to_keep) > 0) {
      this_sel <- which(to_keep == 1)
      this_new <- dplyr::slice(this_prof, this_sel)
      r_new <- nrow(this_new) 
      new_rec[rec_pointer:(rec_pointer + r_new - 1), ] <- this_new
      rec_pointer <- rec_pointer + r_new
    }
  }

  new_rec <- dplyr::slice(new_rec, 1:(rec_pointer - 1)) 
  new_rec[, 2] <- times[new_rec[, 2]]
  
  results <- new("list")
  results$new_rec <- new_rec

  return(results)
}
