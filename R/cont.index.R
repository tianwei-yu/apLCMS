# Computes unique groups
compute_uniq_grp <- function(profile, min.count.run, min.pres = 0.6) {
  grps <- profile  
  ttt <- table(grps)
  ttt <- ttt[ttt >= max(min.count.run * min.pres, 2)]
  unique.grp <- as.numeric(names(ttt))
  return(unique.grp)
}
#' @description
#' Computes the smoothed retention times by using The Nadaraya-Watson kernel regression estimate function.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param times. Retention times vector.
#' @return predicted rt.
#' @examples
#' predict_smoothed_rt(min.run = min.run, times)
predict_smoothed_rt <- function(min.run = 5, times) {
  # ksmooth(x, y, kernel, bandwidth, range, n.points, x.points)
  smooth <- ksmooth(
    seq(-min.run + 1, length(times) + min.run), 
    c(rep(0, min.run),   
    times,       
    rep(0, min.run)),    
    kernel = "box",      
    bandwidth = min.run, 
    x.points = 1:length(times) 
  ) 
  # vector of smoothed estimates for the regression at the corresponding x
  smooth <- smooth$y  
  return(smooth)
}

#' Continuity index
#'
#' This is an internal function. It uses continuity index (or "run filter") to select putative peaks from EIC.
#'
#' @param newprof The matrix containing m/z, retention time, intensity, and EIC label as columns.
#' @param min.pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped
#'  by m/z to be considered a peak.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @return A list is returned.
#' \itemize{
#'   \item new.rec - The matrix containing m/z, retention time, intensity, and EIC label as columns after applying the run filter.
#'   \item height.rec - The vector of peak heights.
#'   \item time.range.rec - The vector of peak retention time span.
#'   \item mz.pres.rec - The vector of proportion of non-missing m/z.
#' }
#' @export
#' @examples
#' cont.index(newprof, min.pres = min.pres, min.run = min.run)
cont.index <- function(newprof,
                       min.pres = 0.6,
                       min.run = 5) {
  
  
  # newprof[,1] : mz
  # newprof[,2] : rt
  # newprof[,3] : intensi
  # newprof[,4] : EIC labels (Extracted Ions Chromatograms)

  # ordering retention time values
  labels <- newprof[,2]
  times <- unique(labels)
  times <- times[order(times)]  
  #time.points <- length(times)
 
  for (i in 1:length(times)) labels[which(newprof[,2] == times[i])] <- i #now labels is the index of time points
  newprof[,2] <- labels  
  
  # set lower bounds of elution time
  min.count.run <- min.run * length(times) / (max(times) - min(times))
  min.run <- round(min.count.run)

  # computes unique groups
  uniq.grp <- compute_uniq_grp(newprof[, 4], min.count.run)

  # rearrange profile data
  newprof <- newprof[newprof[, 4] %in% uniq.grp, ] 
  newprof <- newprof[order(newprof[, 4], newprof[, 1]), ] 

  # computes break points
  breaks <- compute_breaks_3(newprof[, 4])

  # init counters for loop
  new.rec <- newprof * 0
  rec.pointer <- 1
  curr.label <- 1
  height.rec <- mz.pres.rec <- time.range.rec <- rep(0, length(breaks))
  timeline <- rep(0, length(times))

  for (m in 2:length(breaks))
  {
    this.prof <- newprof[(breaks[m - 1] + 1):breaks[m], ]
    this.prof <- this.prof[order(this.prof[, 2]), ]
    
    # init loop variables
    this.mass <- this.prof[, 1]
    this.times <- this.prof[, 2]
    this.intensi <- this.prof[, 3]
    this.grp <- this.prof[1, 4]
    this.timeline <- timeline
    this.timeline[this.times] <- 1
    to.keep <- this.times * 0
    
    # Computes predicted smoothed-rt
    this.smooth <- predict_smoothed_rt(min.run, this.timeline)

    # perform filtering based on the kernel regression estimate
    if (max(this.smooth) >= min.pres) {
      measured.points <- good.points <- timeline
      measured.points[this.times] <- 1

      good.sel <- which(this.smooth >= min.pres)
      good.points[good.sel] <- 1
      for (j in (-min.run):min.run)
      {
        curr.sel <- good.sel + j
        curr.sel <- curr.sel[curr.sel > 0 & curr.sel <= length(times)]
        good.points[curr.sel] <- 1
      }

      measured.points <- measured.points * good.points
      to.keep[which(this.times %in% which(measured.points == 1))] <- 1
    }
    
    # operation over filtered indices 
    if (sum(to.keep) > 0) {
      this.sel <- which(to.keep == 1)
      this.new <- cbind(this.mass[this.sel], this.times[this.sel], this.intensi[this.sel], rep(this.grp, length(this.sel)))
      r.new <- nrow(this.new)
      height.rec[curr.label] <- max(this.intensi[this.sel])
      time.range.rec[curr.label] <- times[max(this.times[this.sel])] - times[min(this.times[this.sel])]
      mz.pres.rec[curr.label] <- length(this.sel) / (max(this.times[this.sel]) - min(this.times[this.sel]) + 1)
      curr.label <- curr.label + 1

      new.rec[rec.pointer:(rec.pointer + r.new - 1), ] <- this.new
      rec.pointer <- rec.pointer + r.new
    }
  }

  new.rec <- new.rec[1:(rec.pointer - 1), ]
  new.rec[, 2] <- times[new.rec[, 2]]

  results <- new("list")
  results$new.rec <- new.rec
  results$height.rec <- height.rec[1:(curr.label - 1)]
  results$time.range.rec <- time.range.rec[1:(curr.label - 1)]
  results$mz.pres.rec <- mz.pres.rec[1:(curr.label - 1)]

  return(results)
}
