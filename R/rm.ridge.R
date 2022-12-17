#' @description
#' Computes retention time intervals of a selected retention time range.
#' @param rt full retention time vector.
#' @param rt_sel Selected retention time range vector.
#' @return A list object:
#' \itemize{
#'   \item over_rt - upper retention time interval
#'   \item under_rt - lower retention time interval
#'   \item within_rt - intermediate retention time interval
#' }
#' @export
compute_rt_intervals <- function(rt, rt_sel){
    rt_max <- max(rt[rt_sel])
    rt_min <- min(rt[rt_sel])

    over_rt <- which(rt > rt_max)
    under_rt <-  which(rt < rt_min)
    within_rt <- which(between(rt, rt_min, rt_max))

    list(over_rt = over_rt, under_rt = under_rt, within_rt = within_rt)
}

#' Removing long ridges at the same m/z.
#' 
#' @description
#' This is an internal function. It substracts a background when an EIC continuously 
#' span more than half the retention time range. The background is estimated through kernel smoothing.
#' @param rt Retention time vector.
#' @param intensity Intensity vector.
#' @param bw Bandwidth for the kernel smoother. A very wide one is used here.
#' @return A vector of intensity values at each rt intervals is returned.
#' @importFrom dplyr between
#' @export
rm.ridge <- function(rt, intensity, bw) {
    this_rt <- which(intensity < quantile(intensity, 0.75))

    rt_intervals <- compute_rt_intervals(rt, this_rt)

    rt_over <- rt_intervals$over_rt 
    rt_under <- rt_intervals$under_rt
    rt_within <- rt_intervals$within_rt

    this.s <- ksmooth(rt[this_rt], intensity[this_rt], x.points = rt[rt_within], kernel = "normal", bandwidth = bw)
    if(sum(is.na(this.s$y)) > 0) return(intensity)
    
    intensity[rt_within] <- intensity[rt_within] - this.s$y
    intensity[rt_over] <- intensity[rt_over] - this.s$y[which(this.s$x == max(this.s$x))[1]]
    intensity[rt_under] <- intensity[rt_under] - this.s$y[which(this.s$x == min(this.s$x))[1]]
    
    intensity[intensity < 0] <- 0
    return(intensity)
}
