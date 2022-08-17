#' @import tibble dplyr
NULL
#> NULL

#' @export
compute_densities <- function(masses, tol, weighted, intensities, bw_func, n = 512) {
  bandwidth <- 0.5 * tol * bw_func(masses)
  if (weighted) {
    weights <- intensities / sum(intensities)
    all.mass.den <- density(masses, weights = weights, bw = bandwidth, n = n)
  } else {
    all.mass.den <- density(masses, bw = bandwidth, n = n)
  }
  return(all.mass.den)
}

#' @export
compute_mass_values <- function(tol, masses, intensi, weighted) {
  n <- 2^min(15, floor(log2(length(masses))) - 2)

  all.mass.den <- compute_densities(masses, tol, weighted, intensi, max, n)

  all.mass.turns <- find.turn.point(all.mass.den$y)
  all.mass.vlys <- all.mass.den$x[all.mass.turns$vlys]
  return(all.mass.vlys)
}

#' @export
compute_breaks <- function(tol, masses, intensi, weighted) {
  all.mass.vlys <- compute_mass_values(tol, masses, intensi, weighted)
  breaks <- c(0, unique(round(approx(masses, 1:length(masses), xout = all.mass.vlys, rule = 2, ties = "ordered")$y))[-1])
  return(breaks)
}


#' @export
increment_counter <- function(pointers, that.n){
  pointers$prof.pointer <- pointers$prof.pointer + that.n
  pointers$height.pointer <- pointers$height.pointer + 1
  pointers$curr.label <- pointers$curr.label + 1

  return(pointers)
}

#' Adaptive binning
#' 
#' This is an internal function. It creates EICs using adaptive binning procedure
#' 
#' @param x A matrix with columns of m/z, retention time, intensity.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be 
#'  considered a peak.
#' @param min.pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped 
#'  by m/z to be considered a peak.
#' @param tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. 
#'  This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy 
#'  level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline.correct After grouping the observations, the highest intensity in each group is found. If the highest 
#'  is lower than this value, the entire group will be deleted. The default value is NA, in which case the program uses the 
#'  75th percentile of the height of the noise groups.
#' @param weighted Whether to weight the local density by signal intensities.
#' @return A list is returned.
#' \itemize{
#'   \item height.rec - The records of the height of each EIC.
#'   \item masses - The vector of m/z values after binning.
#'   \item labels - The vector of retention time after binning.
#'   \item intensi - The vector of intensity values after binning.
#'   \item grps - The EIC labels, i.e. which EIC each observed data point belongs to.
#'   \item times - All the unique retention time values, ordered.
#'   \item tol - The m/z tolerance level.
#'   \item min.count.run - The minimum number of elution time points for a series of signals grouped by m/z to be considered a peak.
#' }
#' @export
#' @examples
#' adaptive.bin(raw.data, min.run = min.run, min.pres = min.pres, tol = tol, baseline.correct = baseline.correct, weighted = intensity.weighted)
adaptive.bin <- function(x,
                         min.run,
                         min.pres,
                         tol,
                         baseline.correct,
                         weighted = FALSE) {
  # order inputs after mz values
  data_table <- tibble::tibble(mz = x$masses, labels = x$labels, intensities = x$intensi) |> dplyr::arrange_at("mz")


  cat(c("m/z tolerance is: ", tol, "\n"))

  times <- x$times #sort(unique(data_table$labels))

  rm(x)
  min_time <- min(times)
  max_time <- max(times)
  time_range <- (max_time - min_time)

  # calculate function parameters
  min.count.run <- min.run * length(times) / time_range
  aver.time.range <- (time_range) / length(times)

  # init data
  newprof <- matrix(0, nrow = length(data_table$mz), ncol = 4)
  height.rec <- matrix(0, nrow = length(data_table$mz), ncol = 3)

  # init counters
  pointers <- list(curr.label = 1, prof.pointer = 1, height.pointer = 1)

  breaks <- compute_breaks(tol, data_table$mz, data_table$intensities, weighted)

  for (i in 1:(length(breaks) - 1))
  {
    
    # get number of scans in bin
    start <- breaks[i] + 1
    end <- breaks[i + 1]

    this_table <- data.frame(labels = data_table$labels[start:end], mz = data_table$mz[start:end], intensities = data_table$intensities[start:end])

    if (length(unique(this_table$labels)) >= min.count.run * min.pres) {
      # reorder in order of labels (scan number)
      this_table <- this_table |> dplyr::arrange_at("labels")
      mass.den <- compute_densities(this_table$mz, tol, weighted, this_table$intensities, median)

      mass.den$y[mass.den$y < min(this_table$intensities) / 10] <- 0
      mass.turns <- find.turn.point(mass.den$y)
      mass.pks <- mass.den$x[mass.turns$pks]
      mass.vlys <- c(-Inf, mass.den$x[mass.turns$vlys], Inf)


      for (j in 1:length(mass.pks))
      {
        # compute boundaries
        boundaries <- compute_boundaries(mass.vlys, mass.pks[j])

        if (length(mass.pks) == 1){
          boundaries$lower <- boundaries$lower - 1
        }

        # get rows which fulfill condition
        that <- this_table |> dplyr::filter(mz > boundaries$lower & mz <= boundaries$upper)

        if (nrow(that) > 0) {
          that <- combine.seq.3(that) |> dplyr::arrange_at("mz")
          that.range <- diff(range(that$labels))

          if (that.range > 0.5 * time_range & length(that$labels) > that.range * min.pres & length(that$labels) / (that.range / aver.time.range) > min.pres) {
            that$intensities <- rm.ridge(that$labels, that$intensities, bw = max(10 * min.run, that.range / 2))

            that <- that |> dplyr::filter(intensities != 0)
          }

          that.n <- length(that$mz)

          newprof[pointers$prof.pointer:(pointers$prof.pointer + that.n - 1), ] <- cbind(that$mz, that$labels, that$intensities, rep(pointers$curr.label, that.n))
          height.rec[pointers$height.pointer, ] <- c(pointers$curr.label, that.n, max(that$intensities))

          # increment counters
          pointers <- increment_counter(pointers, that.n)
        }
      }
    } else {
      if (runif(1) < 0.05) {
        this_table <- this_table |> dplyr::arrange_at("labels")

        that.merged <- combine.seq.3(this_table)
        that.n <- nrow(that.merged)

        newprof[pointers$prof.pointer:(pointers$prof.pointer + that.n - 1), ] <- cbind(that.merged$mz, that.merged$labels, that.merged$intensities, rep(pointers$curr.label, that.n))
        height.rec[pointers$height.pointer, ] <- c(pointers$curr.label, that.n, max(that.merged$intensities))

        # increment counters
        pointers <- increment_counter(pointers, that.n)
      }
    }
  }

  newprof <- newprof[1:(pointers$prof.pointer - 1), ]
  height.rec <- height.rec[1:(pointers$height.pointer - 1), ]

  newprof <- newprof[order(newprof[, 1], newprof[, 2]), ]

  raw.prof <- new("list")
  raw.prof$height.rec <- height.rec
  raw.prof$masses <- newprof[, 1]
  raw.prof$labels <- newprof[, 2]
  raw.prof$intensi <- newprof[, 3]
  raw.prof$grps <- newprof[, 4]
  raw.prof$times <- times
  raw.prof$tol <- tol
  raw.prof$min.count.run <- min.count.run

  return(raw.prof)
}
