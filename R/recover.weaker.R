#' @export
duplicate.row.remove <- function(new.table) {
  new.table <- new.table[order(new.table[, 1], new.table[, 2], new.table[, 5]), ]
  n <- 1
  m <- 2
  to.remove <- rep(0, nrow(new.table))

  while (m <= nrow(new.table)) {
    if (abs(new.table[m, 1] - new.table[n, 1]) < 1e-10 &
      abs(new.table[m, 2] - new.table[n, 2]) < 1e-10 &
      abs(new.table[m, 5] - new.table[n, 5]) < 1e-10) {
      to.remove[m] <- 1
      m <- m + 1
    } else {
      n <- m
      m <- m + 1
    }
    # cat("*(", n, m, ")")
  }

  if (sum(to.remove) > 0) new.table <- new.table[-which(to.remove == 1), ]
  new.table
}

#' Compute all time values from base curve.
#' @param base_curve Basis curve
#' @export
compute_all_times <- function(base_curve) {
  all_times <- base_curve[, 1]
  if (all_times[1] > 0) all_times <- c(0, all_times)
  all_times <- c(all_times, 2 * all_times[length(all_times)] - all_times[length(all_times) - 1])
  all_times <- (all_times[1:(length(all_times) - 1)] + all_times[2:length(all_times)]) / 2
  all_times <- all_times[2:length(all_times)] - all_times[1:(length(all_times) - 1)]
  return(all_times)
}

#' @export
compute_base_curve <- function(x) {
  base_curve <- unique(x)
  base_curve <- base_curve[order(base_curve)]
  base_curve <- cbind(base_curve, base_curve * 0)
  return(base_curve)
}


#' Normalize vector so that sum(vec) = 1
l2normalize <- function(x) {
  x / sum(x)
}

#' @export
compute_mass_density <- function(mz,
                                 intensities,
                                 bandwidth,
                                 intensity_weighted) {
  if (intensity_weighted) {
    mass_density <- density(
      mz,
      weights = l2normalize(intensities),
      bw = bandwidth
    )
  } else {
    mass_density <- density(mz, bw = bandwidth)
  }
  return(mass_density)
}

#' @export
get_custom_chr_tol <- function(use.observed.range,
                               pk.times,
                               chr.range,
                               aligned.ftrs) {
  custom.chr.tol <- rep(chr.range, nrow(aligned.ftrs))

  if (use.observed.range) {
    # check observed rt range across ALL SAMPLES
    all_peak_rts <- pk.times[, 5:ncol(pk.times)]
    observed.chr.range <- (apply(all_peak_rts, 1, max) - apply(all_peak_rts, 1, min)) / 2
    sufficient_rts <- apply(!is.na(all_peak_rts), 1, sum) >= 5
    selection <- which(sufficient_rts & custom.chr.tol > observed.chr.range)
    custom.chr.tol[selection] <- observed.chr.range[selection]
  }

  return(custom.chr.tol)
}

#' @export
compute_target_time <- function(aligned_rts, orig.time, adjusted.time) {
  to.use <- get_times_to_use(orig.time, adjusted.time)
  orig.time <- orig.time[to.use]
  adjusted.time <- adjusted.time[to.use]

  sel.non.na <- which(!is.na(aligned_rts))
  if (length(adjusted.time) >= 4) {
    sp <- interpSpline(orig.time ~ adjusted.time, na.action = na.omit)
    aligned_rts[sel.non.na] <- predict(sp, aligned_rts[sel.non.na])$y
  }
}

#' @export
get_times_to_use <- function(orig.time, adjusted.time) {
  ttt.0 <- table(orig.time)
  ttt <- table(adjusted.time)
  to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt == 1]) & orig.time %in% as.numeric(names(ttt.0)[ttt.0 == 1]))
  if (length(to.use) > 2000) {
    to.use <- sample(to.use, 2000, replace = FALSE)
  }
  return(to.use)
}

#' @export
compute_breaks_2 <- function(data_table, orig.tol) {
  all.mass.den <- density(
    data_table$mz,
    weights = l2normalize(data_table$intensities),
    bw = 0.5 * orig.tol * max(data_table$mz),
    n = 2^min(15, floor(log2(length(data_table$mz))) - 2)
  )

  all.mass.turns <- find.turn.point(all.mass.den$y)
  all.mass.vlys <- all.mass.den$x[all.mass.turns$vlys]
  breaks <- c(0, unique(round(approx(data_table$mz, seq_along(data_table$mz), xout = all.mass.vlys, rule = 2, ties = "ordered")$y))[-1])

  return(breaks)
}

#' @export
get_mzrange_bound_indices <- function(aligned.ftrs, masses, breaks, i, custom.mz.tol) {
  if (aligned.ftrs[i, "mz"] <= masses[breaks[2]]) {
    this.found <- c(1, 2)
  } else {
    this.found <- c(which(abs(masses[breaks] - aligned.ftrs[i, "mz"]) < custom.mz.tol[i]), min(which(masses[breaks] > aligned.ftrs[i, "mz"])), max(which(masses[breaks] < aligned.ftrs[i, "mz"]))) + 1
    this.found <- c(min(this.found), max(this.found))
  }
  return(this.found)
}

#' @export
get_raw_features_in_mzrange <- function(data_table, aligned.ftrs, breaks, i, custom.mz.tol) {
  this.found <- get_mzrange_bound_indices(aligned.ftrs, data_table$mz, breaks, i, custom.mz.tol)
  this.sel <- (breaks[this.found[1]] + 1):breaks[this.found[2]]
  features <- data_table |> dplyr::slice(this.sel)
  return(features)
}

#' Recover weak signals in some profiles that is not identified as a peak, but corresponds to identified peaks in other spectra.
#' 
#' @description
#' Given the aligned feature table, some features are identified in a subgroup of spectra. This doesn't mean they don't exist in the other spectra. 
#' The signal could be too low to pass the run filter. Thus after obtaining the aligned feature table, this function re-analyzes each spectrum to 
#' try and fill in the holes in the aligned feature table.
#' @param filename the cdf file name from which weaker signal is to be recovered.
#' @param loc the location of the filename in the vector of filenames.
#' @param aligned.ftrs matrix, with columns of m/z values, elution times, signal strengths in each spectrum.
#' @param pk.times matrix, with columns of m/z, median elution time, and elution times in each spectrum.
#' @param align.mz.tol the m/z tolerance used in the alignment.
#' @param align.chr.tol the elution time tolerance in the alignment.
#' @param this.f1 The matrix which is the output from proc.to.feature().
#' @param this.f2 The matrix which is the output from proc.to.feature(). The retention time in this object have been adjusted by the function adjust.time().
#' @param mz.range The m/z around the feature m/z to search for observations. The default value is NA, in which case 1.5 times the m/z tolerance in the aligned object will be used.
#' @param chr.range The retention time around the feature retention time to search for observations. The default value is NA, in which case 0.5 times the retention time tolerance in the aligned object will be used.
#' @param use.observed.range If the value is TRUE, the actual range of the observed locations of the feature in all the spectra will be used.
#' @param orig.tol The mz.tol parameter provided to the proc.cdf() function. This helps retrieve the intermediate file.
#' @param min.bw The minimum bandwidth to use in the kernel smoother.
#' @param max.bw The maximum bandwidth to use in the kernel smoother.
#' @param bandwidth A value between zero and one. Multiplying this value to the length of the signal along the time axis helps determine the bandwidth in the kernel smoother used for peak identification.
#' @param recover.min.count minimum number of raw data points to support a recovery.
#' @param intensity.weighted Whether to use intensity to weight mass density estimation.
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item aligned.ftrs - A matrix, with columns of m/z values, elution times, and signal strengths in each spectrum.
#'   \item pk.times - A matrix, with columns of m/z, median elution time, and elution times in each spectrum.
#'   \item mz.tol - The m/z tolerance in the aligned object.
#'   \item chr.tol - The elution time tolerance in the aligned object.
#' }
#' @export
#' @examples
#' recover.weaker(filename, loc, aligned.ftrs, pk.times, align.mz.tol, align.chr.tol, this.f1, this.f2)
recover.weaker <- function(filename,
                           sample_name,
                           aligned.ftrs,
                           pk.times,
                           align.mz.tol,
                           align.chr.tol,
                           this.f1,
                           this.f2,
                           mz.range = NA,
                           chr.range = NA,
                           use.observed.range = TRUE,
                           orig.tol = 1e-5,
                           min.bw = NA,
                           max.bw = NA,
                           bandwidth = .5,
                           recover.min.count = 3,
                           intensity.weighted = FALSE) {

  # load raw data
  this.raw <- load_file(filename)
  times <- this.raw$times
  data_table <- tibble::tibble(mz = this.raw$masses, labels = this.raw$labels, intensities = this.raw$intensi) |> dplyr::arrange_at("mz")
  rm(this.raw)
  masses <- data_table$mz

  # Initialize parameters with default values
  if (is.na(mz.range)) mz.range <- 1.5 * align.mz.tol
  if (is.na(chr.range)) chr.range <- align.chr.tol / 2
  if (is.na(min.bw)) min.bw <- diff(range(times, na.rm = TRUE)) / 60
  if (is.na(max.bw)) max.bw <- diff(range(times, na.rm = TRUE)) / 15
  if (min.bw >= max.bw) min.bw <- max.bw / 4


  base.curve <- compute_base_curve(sort(times))
  aver.diff <- mean(diff(base.curve[, 1]))
  all.times <- compute_all_times(base.curve)

  this.ftrs <- aligned.ftrs[, sample_name]
  this.times <- pk.times[, sample_name]

  custom.mz.tol <- mz.range * aligned.ftrs$mz
  custom.chr.tol <- get_custom_chr_tol(
    use.observed.range,
    pk.times,
    chr.range,
    aligned.ftrs
  )

  target.time <- compute_target_time(
    aligned.ftrs[, "rt"],
    round(this.f1[, "pos"], 5),
    round(this.f2[, "pos"], 5)
  )

  breaks <- compute_breaks_2(data_table, orig.tol)

  this.mz <- rep(NA, length(this.ftrs))

  for (i in seq_along(this.ftrs))
  {
    if (this.ftrs[i] == 0 && aligned.ftrs[i, "mz"] < max(masses)) {
      features <- get_raw_features_in_mzrange(data_table, aligned.ftrs, breaks, i, custom.mz.tol)

      mass.den <- compute_mass_density(
        mz = features$mz,
        intensities = features$intensities,
        bandwidth = 0.5 * orig.tol * aligned.ftrs[i, "mz"],
        intensity_weighted = intensity.weighted
      )

      # find peaks in mz range in raw data
      mass.turns <- find.turn.point(mass.den$y)
      mass.pks <- mass.den$x[mass.turns$pks] # mz values with highest density
      mass.vlys <- c(-Inf, mass.den$x[mass.turns$vlys], Inf) # masses with lowest densities values -> valley
      mass.pks <- mass.pks[which(abs(mass.pks - aligned.ftrs[i, "mz"]) < custom.mz.tol[i] / 1.5)]

      if (length(mass.pks) > 0) {
        this.rec <- matrix(c(Inf, Inf, Inf), nrow = 1)
        for (peak in mass.pks)
        {
          # get mass values of valleys the closest to the peak
          mass.lower <- max(mass.vlys[mass.vlys < peak])
          mass.upper <- min(mass.vlys[mass.vlys > peak])

          that <- features |>
            dplyr::filter(mz > mass.lower && mz <= mass.upper) |>
            dplyr::arrange_at("labels")

          if (nrow(that) > recover.min.count) {
            that.prof <- combine.seq.3(that$labels, that$mz, that$intensities)
            that.mass <- sum(that.prof[, 1] * that.prof[, 3]) / sum(that.prof[, 3])
            curr.rec <- c(that.mass, NA, NA)

            if (nrow(that.prof) < 10) {
              if (!is.na(target.time[i])) {
                thee.sel <- which(abs(that.prof[, 2] - target.time[i]) < custom.chr.tol[i] * 2)
              } else {
                thee.sel <- 1:nrow(that.prof)
              }

              if (length(thee.sel) > recover.min.count) {
                if (length(thee.sel) > 1) {
                  that.inte <- interpol.area(that.prof[thee.sel, 2], that.prof[thee.sel, 3], base.curve[, 1], all.times)
                } else {
                  that.inte <- that.prof[thee.sel, 3] * aver.diff
                }
                curr.rec[3] <- that.inte
                curr.rec[2] <- median(that.prof[thee.sel, 2])
                this.rec <- rbind(this.rec, curr.rec)
              }
            } else {
              this <- that.prof[, 2:3]
              this <- this[order(this[, 1]), ]
              this.span <- range(this[, 1])
              this.curve <- base.curve[base.curve[, 1] >= this.span[1] & base.curve[, 1] <= this.span[2], ]
              this.curve[this.curve[, 1] %in% this[, 1], 2] <- this[, 2]

              bw <- min(max(bandwidth * (max(this[, 1]) - min(this[, 1])), min.bw), max.bw)
              this.smooth <- ksmooth(this.curve[, 1], this.curve[, 2], kernel = "normal", bandwidth = bw)
              smooth.y <- this.smooth$y
              turns <- find.turn.point(smooth.y)
              pks <- this.smooth$x[turns$pks]
              vlys <- this.smooth$x[turns$vlys]
              vlys <- c(-Inf, vlys, Inf)

              pks.n <- pks
              for (m in 1:length(pks))
              {
                this.vlys <- c(max(vlys[which(vlys < pks[m])]), min(vlys[which(vlys > pks[m])]))
                pks.n[m] <- sum(this[, 1] >= this.vlys[1] & this[, 1] <= this.vlys[2])
              }

              if (!is.na(target.time[i])) {
                pks.d <- abs(pks - target.time[i]) # distance from the target peak location
                pks.d[pks.n == 0] <- Inf
                pks <- pks[which(pks.d == min(pks.d))[1]]
              } else {
                pks <- pks[pks.n > recover.min.count]
              }

              all.pks <- pks
              all.vlys <- vlys
              all.this <- this

              if (length(all.pks) > 0) {
                for (pks.i in 1:length(all.pks))
                {
                  pks <- all.pks[pks.i]
                  vlys <- c(max(all.vlys[which(all.vlys < pks)]), min(all.vlys[which(all.vlys > pks)]))

                  this <- all.this[which(all.this[, 1] >= vlys[1] & all.this[, 1] <= vlys[2]), ]
                  if (is.null(nrow(this))) {
                    curr.rec[3] <- this[2] * aver.diff
                    curr.rec[2] <- this[1]
                  } else {
                    x <- this[, 1]
                    y <- this[, 2]

                    if (nrow(this) >= 10) {
                      miu <- sum(x * y) / sum(y)
                      sigma <- sqrt(sum(y * (x - miu)^2) / sum(y))
                      if (sigma == 0) {
                        curr.rec[3] <- sum(y) * aver.diff
                        curr.rec[2] <- miu
                      } else {
                        fitted <- dnorm(x, mean = miu, sd = sigma)
                        this.sel <- y > 0 & fitted / dnorm(miu, mean = miu, sd = sigma) > 1e-2
                        sc <- exp(sum(fitted[this.sel]^2 * log(y[this.sel] / fitted[this.sel]) / sum(fitted[this.sel]^2)))
                      }
                    } else {
                      sc <- interpol.area(x, y, base.curve[, 1], all.times)
                      miu <- median(x)
                    }
                    curr.rec[3] <- sc
                    curr.rec[2] <- miu
                  }
                  this.rec <- rbind(this.rec, curr.rec)
                }
              }
            }
          }
        }

        if (!is.na(target.time[i])) {
          this.sel <- which(abs(this.rec[, 2] - target.time[i]) < custom.chr.tol[i])
        } else {
          this.sel <- 1:nrow(this.rec)
          this.sel <- this.sel[this.sel != 1]
        }


        if (length(this.sel) > 0) {
          if (length(this.sel) > 1) {
            if (!is.na(target.time[i])) {
              this.d <- (this.rec[, 2] - target.time[i])^2 / custom.chr.tol[i]^2 + (this.rec[, 1] - aligned.ftrs[i, 1])^2 / custom.mz.tol[i]^2
              this.sel <- which(this.d == min(this.d))[1]
            } else {
              this.d <- abs(this.rec[, 1] - aligned.ftrs[i, 1])
              this.sel <- which(this.d == min(this.d))[1]
            }
          }
          this.pos.diff <- abs(this.f1[, 2] - this.rec[this.sel, 2])
          this.pos.diff <- which(this.pos.diff == min(this.pos.diff))[1]
          this.f1 <- rbind(this.f1, c(this.rec[this.sel, 1], this.rec[this.sel, 2], NA, NA, this.rec[this.sel, 3]))
          this.time.adjust <- (-this.f1[this.pos.diff, 2] + this.f2[this.pos.diff, 2])
          this.f2 <- rbind(
            this.f2,
            c(
              this.rec[this.sel, 1],
              this.rec[this.sel, 2] + this.time.adjust,
              NA,
              NA,
              this.rec[this.sel, 3],
              grep(sample_name, colnames(aligned.ftrs)),
              NA
            )
          )
          this.ftrs[i] <- this.rec[this.sel, 3]
          this.times[i] <- this.rec[this.sel, 2] + this.time.adjust
          this.mz[i] <- this.rec[this.sel, 1]
        }
      }
    }
  }
  to.return <- new("list")
  to.return$this.mz <- this.mz
  to.return$this.ftrs <- this.ftrs
  to.return$this.times <- this.times
  to.return$this.f1 <- duplicate.row.remove(this.f1)
  to.return$this.f2 <- duplicate.row.remove(this.f2)

  return(to.return)
}
