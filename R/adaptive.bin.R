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
#' @examples
#' adaptive.bin(raw.data, min.run = min.run, min.pres = min.pres, tol = tol, baseline.correct = baseline.correct, weighted = intensity.weighted)
adaptive.bin <- function(x,
                         min.run,
                         min.pres,
                         tol,
                         baseline.correct,
                         weighted = FALSE) {
  # order inputs after mz values
  masses <- x$masses
  labels <- x$labels
  intensi <- x$intensi
  times <- x$times

  rm(x)

  curr.order <- order(masses)
  intensi <- intensi[curr.order]
  labels <- labels[curr.order]
  masses <- masses[curr.order]

  rm(curr.order)

  cat(c("m/z tolerance is: ", tol, "\n"))

  times <- unique(labels)
  times <- times[order(times)]

  # calculate function parameters
  min.count.run <- min.run * length(times) / (max(times) - min(times))
  time.range <- diff(range(times))
  aver.time.range <- (max(labels) - min(labels)) / length(times)

  # init data
  newprof <- matrix(0, nrow = length(masses), ncol = 4)
  height.rec <- matrix(0, nrow = length(masses), ncol = 3)

  # init counters
  curr.label <- 1
  prof.pointer <- 1
  height.pointer <- 1

  breaks <- compute_breaks(tol, masses, intensi, weighted)

  for (i in 1:(length(breaks) - 1))
  {
    start <- breaks[i] + 1
    end <- breaks[i + 1]
    # get number of scans in bin
    this.labels <- labels[start: end]

    if (length(unique(this.labels)) >= min.count.run * min.pres) {
      # extract mz and intensity values for bin
      this.masses <- masses[start:end]
      this.intensi <- intensi[start:end]

      # reorder in order of labels (scan number)
      curr.order <- order(this.labels)
      this.masses <- this.masses[curr.order]
      this.intensi <- this.intensi[curr.order]
      this.labels <- this.labels[curr.order]

      mass.den <- compute_densities(this.masses, tol, weighted, this.intensi, median)

      mass.den$y[mass.den$y < min(this.intensi) / 10] <- 0
      mass.turns <- find.turn.point(mass.den$y)
      mass.pks <- mass.den$x[mass.turns$pks]
      mass.vlys <- c(-Inf, mass.den$x[mass.turns$vlys], Inf)


      for (j in 1:length(mass.pks))
      {
        # compute boundaries
        mass.lower <- max(mass.vlys[mass.vlys < mass.pks[j]])
        mass.upper <- min(mass.vlys[mass.vlys > mass.pks[j]])

        if (length(mass.pks) == 1) mass.lower <- mass.lower - 1

        # compute if we are in mass range from mass.lower to mass.upper
        mass.sel <- which(this.masses > mass.lower & this.masses <= mass.upper)

        if (length(mass.sel) > 0) {

          # get rows which fulfill condition
          that.labels <- this.labels[mass.sel]
          that.masses <- this.masses[mass.sel]
          that.intensi <- this.intensi[mass.sel]

          # rearrange in order of labels
          that.merged <- combine.seq.3(that.labels, that.masses, that.intensi)
          if (nrow(that.merged) == 1) {
            new.merged <- that.merged
          } else {
            new.merged <- that.merged[order(that.merged[, 1]), ]
          }

          that.labels <- new.merged[, 2]
          that.masses <- new.merged[, 1]
          that.intensi <- new.merged[, 3]
          that.range <- diff(range(that.labels))

          if (that.range > 0.5 * time.range & length(that.labels) > that.range * min.pres & length(that.labels) / (diff(range(that.labels)) / aver.time.range) > min.pres) {
            that.intensi <- rm.ridge(that.labels, that.intensi, bw = max(10 * min.run, that.range / 2))

            # filter out 0 entries
            that.sel <- which(that.intensi != 0)
            that.labels <- that.labels[that.sel]
            that.masses <- that.masses[that.sel]
            that.intensi <- that.intensi[that.sel]
          }

          that.n <- length(that.masses)

          newprof[prof.pointer:(prof.pointer + that.n - 1), ] <- cbind(that.masses, that.labels, that.intensi, rep(curr.label, that.n))
          height.rec[height.pointer, ] <- c(curr.label, that.n, max(that.intensi))

          # increment counters
          prof.pointer <- prof.pointer + that.n
          height.pointer <- height.pointer + 1
          curr.label <- curr.label + 1
        }
      }
    } else {
      if (runif(1) < 0.05) {

        # reassignment
        this.masses <- masses[start:end]
        this.intensi <- intensi[start:end]

        # reordering
        curr.order <- order(this.labels)
        this.masses <- this.masses[curr.order]
        this.intensi <- this.intensi[curr.order]
        this.labels <- this.labels[curr.order]

        that.merged <- combine.seq.3(this.labels, this.masses, this.intensi)
        that.n <- nrow(that.merged)

        newprof[prof.pointer:(prof.pointer + that.n - 1), ] <- cbind(that.merged, rep(curr.label, that.n))
        height.rec[height.pointer, ] <- c(curr.label, that.n, max(that.merged[, 3]))

        # increment counters
        prof.pointer <- prof.pointer + that.n
        height.pointer <- height.pointer + 1
        curr.label <- curr.label + 1
      }
    }
  }

  newprof <- newprof[1:(prof.pointer - 1), ]
  height.rec <- height.rec[1:(height.pointer - 1), ]

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
