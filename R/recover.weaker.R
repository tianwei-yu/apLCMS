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

compute_mass_density <- function(mz, intensities, bandwidth, intensity_weighted) {
  if (intensity_weighted) {
    mass_density <- density(mz, weights = l2normalize(intensities), bw = bandwidth)
  } else {
    mass_density <- density(mz, bw = bandwidth)
  }
  return(mass_density)
}

recover.weaker <- function(
  filename, loc, aligned.ftrs, pk.times, align.mz.tol, align.chr.tol, this.f1, this.f2, mz.range = NA, chr.range = NA, use.observed.range = TRUE, orig.tol = 1e-5, min.bw = NA, max.bw = NA, bandwidth = .5, recover.min.count = 3, intensity.weighted = FALSE) {
  
  # load raw data
  this.raw <- load_file(filename)
  masses <- this.raw$masses
  intensi <- this.raw$intensi
  labels <- this.raw$labels
  times <- this.raw$times
  rm(this.raw)

  # Initialize parameters with default values
  if (is.na(mz.range)) mz.range <- 1.5 * align.mz.tol
  if (is.na(chr.range)) chr.range <- align.chr.tol / 2
  if (is.na(min.bw)) min.bw <- diff(range(times, na.rm = TRUE)) / 60
  if (is.na(max.bw)) max.bw <- diff(range(times, na.rm = TRUE)) / 15
  if (min.bw >= max.bw) min.bw <- max.bw / 4


  masses <- c(masses, -100000)
  mass.breaks <- which(masses[1:(length(masses) - 1)] > masses[2:length(masses)])
  mass.breaks <- c(0, mass.breaks)
  masses <- masses[1:(length(masses) - 1)]

  data_table <- tibble::tibble(mz = masses, labels = labels, intensities = intensi) |> dplyr::arrange_at("mz")

  masses <- data_table$mz

  base.curve <- compute_base_curve(sort(times))
  aver.diff <- mean(diff(base.curve[,1]))
  all.times <- compute_all_times(base.curve)

  this.ftrs <- aligned.ftrs[, (loc + 4)]
  this.times <- pk.times[, (loc + 4)]

  custom.mz.tol <- mz.range * aligned.ftrs[, 1]
  custom.chr.tol <- rep(chr.range, nrow(aligned.ftrs))

  if (use.observed.range) {
    observed.chr.range <- (apply(pk.times[, 5:ncol(pk.times)], 1, max) - apply(pk.times[, 5:ncol(pk.times)], 1, min)) / 2
    num.present <- apply(!is.na(pk.times[, 5:ncol(pk.times)]), 1, sum)
    custom.chr.tol[which(num.present >= 5 & custom.chr.tol > observed.chr.range)] <- observed.chr.range[which(num.present >= 5 & custom.chr.tol > observed.chr.range)]
  }
  orig.time <- round(this.f1[, 2], 5)
  adjusted.time <- round(this.f2[, 2], 5)
  ttt.0 <- table(orig.time)
  ttt <- table(adjusted.time)
  
  to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt == 1]) & orig.time %in% as.numeric(names(ttt.0)[ttt.0 == 1]))
  if (length(to.use) > 2000) to.use <- sample(to.use, 2000, replace = FALSE)

  sel.non.na <- which(!is.na(aligned.ftrs[, 2]))
  target.time <- aligned.ftrs[, 2]

  orig.time <- orig.time[to.use]
  adjusted.time <- adjusted.time[to.use]

  if (length(adjusted.time) >= 4) {
    sp <- interpSpline(orig.time ~ adjusted.time, na.action = na.omit)
    target.time[sel.non.na] <- predict(sp, aligned.ftrs[sel.non.na, 2])$y
  }

  l <- length(masses)
  curr.bw <- 0.5 * orig.tol * max(masses)
  all.mass.den <- density(
    masses,
    weights = l2normalize(data_table$intensities),
    bw = curr.bw,
    n = 2^min(15, floor(log2(l)) - 2)
  )

  all.mass.turns <- find.turn.point(all.mass.den$y)
  all.mass.vlys <- all.mass.den$x[all.mass.turns$vlys]
  breaks <- c(0, unique(round(approx(masses, 1:l, xout = all.mass.vlys, rule = 2, ties = "ordered")$y))[-1])
  this.mz <- rep(NA, length(this.ftrs))

  for (i in 1:length(this.ftrs))
  {
    if (this.ftrs[i] == 0 & aligned.ftrs[i, 1] < masses[breaks[length(breaks)]]) {
      if (aligned.ftrs[i, 1] <= masses[breaks[2]]) {
        this.found <- c(1, 2)
      } else {
        this.found <- c(which(abs(masses[breaks] - aligned.ftrs[i, 1]) < custom.mz.tol[i]), min(which(masses[breaks] > aligned.ftrs[i, 1])), max(which(masses[breaks] < aligned.ftrs[i, 1]))) + 1
        this.found <- c(min(this.found), max(this.found))
      }

      if (length(this.found) > 1) {
        this.sel <- (breaks[this.found[1]] + 1):breaks[this.found[2]]
        this.masses <- masses[this.sel]
        this.labels <- data_table$labels[this.sel]
        this.intensi <- data_table$intensities[this.sel]

        mass.den <- compute_mass_density(
          mz = this.masses,
          intensities = this.intensi,
          bandwidth = 0.5 * orig.tol * aligned.ftrs[i, 1],
          intensity_weighted = intensity.weighted
        )

        mass.turns <- find.turn.point(mass.den$y)
        mass.pks <- mass.den$x[mass.turns$pks]
        mass.vlys <- c(-Inf, mass.den$x[mass.turns$vlys], Inf)
        mass.pks <- mass.pks[which(abs(mass.pks - aligned.ftrs[i, 1]) < custom.mz.tol[i] / 1.5)]

        if (length(mass.pks) > 0) {
          this.rec <- matrix(c(Inf, Inf, Inf), nrow = 1)
          for (k in 1:length(mass.pks))
          {
            mass.lower <- max(mass.vlys[mass.vlys < mass.pks[k]])
            mass.upper <- min(mass.vlys[mass.vlys > mass.pks[k]])

            that.sel <- which(this.masses > mass.lower & this.masses <= mass.upper)
            if (length(that.sel) > recover.min.count) {
              that.labels <- this.labels[that.sel]
              that.masses <- this.masses[that.sel]
              that.intensi <- this.intensi[that.sel]

              that.order <- order(that.labels)
              that.labels <- that.labels[that.order]
              that.masses <- that.masses[that.order]
              that.intensi <- that.intensi[that.order]

              that.prof <- combine.seq.3(that.labels, that.masses, that.intensi)

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
            this.f2 <- rbind(this.f2, c(this.rec[this.sel, 1], this.rec[this.sel, 2] + this.time.adjust, NA, NA, this.rec[this.sel, 3], loc, NA))
            this.ftrs[i] <- this.rec[this.sel, 3]
            this.times[i] <- this.rec[this.sel, 2] + this.time.adjust
            this.mz[i] <- this.rec[this.sel, 1]
          }
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
