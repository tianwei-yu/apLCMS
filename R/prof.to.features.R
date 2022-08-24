
#' Validate that provided inputs match expected, exit execution otherwise.
#' @param shape.model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian".
#'  When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param estim.method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @export
validate_inputs <- function(shape.model, estim.method) {
  if (!shape.model %in% c("Gaussian", "bi-Gaussian")) {
    stop("shape.model argument must be 'Gaussian' or 'bi-Gaussian'")
  }
  if (!estim.method %in% c("moment", "EM")) {
    stop("estim.method argument must be 'moment' or 'EM'")
  }
}

#' Initialize minimum and maximum bandwidth values if none given. Ensure that minimum bandwidth is lower that maximum, else set minimum to 1/4 of maximum value.
#' @param min.bw The minimum bandwidth to use in the kernel smoother.
#' @param max.bw The maximum bandwidth to use in the kernel smoother.
#' @param profile Profile table with shape number-of-features*4. The table contains following columns:
#' \itemize{
#'   \item mz - float - mass-to-charge ratio of feature
#'   \item rt - float - retention time of features
#'   \item intensity - float - intensity of features
#'   \item group_number - integer - group number assigned to each feature based on their rt similarity
#' }
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item min.bw - float - Minimum bandwidth.
#'   \item max.bw - float - Maximum bandwidth
#' @export
preprocess_bandwidth <- function(min.bw, max.bw, profile) {
  if (is.na(min.bw)) {
    min.bw <- diff(range(profile[, 2], na.rm = TRUE)) / 60
  }
  if (is.na(max.bw)) {
    max.bw <- diff(range(profile[, 2], na.rm = TRUE)) / 15
  }
  if (min.bw >= max.bw) {
    min.bw <- max.bw / 4
  }

  return(list("min.bw" = min.bw, "max.bw" = max.bw))
}

#' Convert input matrix to a dataframe with column names (see source code for the names).
#' @param profile Profile table with shape number-of-features*4. The table contains following columns:
#' \itemize{
#'   \item float - mass-to-charge ratio of feature
#'   \item float - retention time of features
#'   \item float - intensity of features
#'   \item integer - group number assigned to each feature based on their rt similarity
#' }
#' @return  Returns a dataframe with shape number-of-features*4. The columns are as follows:
#' \itemize{
#'   \item mz - float - mass-to-charge ratio of feature
#'   \item rt - float - retention time of features
#'   \item intensity - float - intensity of features
#'   \item group_number - integer - group number assigned to each feature based on their rt similarity
#' }
#' @export
preprocess_profile <- function(profile) {
  keys <- c("mz", "rt", "intensity", "group_number")
  colnames(profile) <- keys

  return(data.frame(profile))
}

#' Compute parameters of chromatographic peak shape if peaks are considered to be gaussian
#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param bw Bandwidth vector to use in the kernel smoother.
#' @param component.eliminate When a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param BIC.factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#'  models with more peaks are penalized more.
#' @param aver_diff Average retention time difference across RTs of all features.
#' @return Returns a single-row vector or a table object with the following items/columns:
#' \itemize{
#'   \item miu - float - mean value of the gaussian curve
#'   \item sigma - float - standard deviation of the gaussian curve
#'   \item sigma - float - standard deviation of the gaussian curve
#'   \item scale - float - estimated total signal strength (total area of the estimated normal curve)
#'}
#' @export
compute_gaussian_peak_shape <- function(rt_profile, power, bw, component.eliminate, BIC.factor, aver_diff) {
  rt_peak_shape <- normix.bic(rt_profile[, "base.curve"], rt_profile[, 2], power = power, bw = bw, eliminate = component.eliminate, BIC.factor = BIC.factor, aver_diff = aver_diff)$param
  if (nrow(rt_peak_shape) == 1) {
    rt_peak_shape <- c(rt_peak_shape[1, 1:2], rt_peak_shape[1, 2], rt_peak_shape[1, 3])
  } else {
    rt_peak_shape <- cbind(rt_peak_shape[, 1:2], rt_peak_shape[, 2], rt_peak_shape[, 3])
  }
  return(rt_peak_shape)
}

#' This function solves the value of a using the x, t, a from the previous step, and sigma.1, and sigma.2 (original authors' comment).
#' @export
solve.a <- function(x, t, a, sigma.1, sigma.2) {
  # This function is a part of bigauss.esti.EM and is not covered by any of test-cases
  w <- x * (as.numeric(t < a) / sigma.1 + as.numeric(t >= a) / sigma.2)
  return(sum(t * w) / sum(w))
}

#' This function prepares the parameters required for latter compuation. u, v, and sum of x (original authors' comment).
#' @export
prep.uv <- function(x, t, a) {
  # This function is a part of bigauss.esti.EM and is not covered by any of test-cases
  temp <- (t - a)^2 * x
  u <- sum(temp * as.numeric(t < a))
  v <- sum(temp * as.numeric(t >= a))
  return(list(
    u = u,
    v = v,
    x.sum = sum(x)
  ))
}

#' This function takes the value intensity level x, retention time t and assumed breaking point a, calculates the square estimated of sigma.1 and sigma.2 (original authors' comment).
#' @export
solve.sigma <- function(x, t, a) {
  # This function is a part of bigauss.esti.EM and is not covered by any of test-cases
  tt <- prep.uv(x, t, a)
  sigma.1 <- tt$u / tt$x.sum * ((tt$v / tt$u)^(1 / 3) + 1)
  sigma.2 <- tt$v / tt$x.sum * ((tt$u / tt$v)^(1 / 3) + 1)
  return(list(
    sigma.1 = sigma.1,
    sigma.2 = sigma.2
  ))
}

#' @description
#' Function takes into x and t, and then computes the value of sigma.1, sigma.2 and a using iterative method. the returned values include estimated sigmas,
#' a and a boolean variable on whether the termination criteria is satified upon the end of the program (original authors' comment).
#' @export
bigauss.esti.EM <- function(t, x, max.iter = 50, epsilon = 0.005, power = 1, do.plot = FALSE, truth = NA, sigma.ratio.lim = c(0.3, 1)) {
  # This function is not covered by any test case
  sel <- which(x > 1e-10)
  if (length(sel) == 0) {
    return(c(median(t), 1, 1, 0))
  }
  if (length(sel) == 1) {
    return(c(t[sel], 1, 1, 0))
  }
  t <- t[sel]
  x <- x[sel]

  ## epsilon is the threshold for continuing the iteration. change in
  ## a smaller than epsilon will terminate the iteration.
  ## epsilon <- min(diff(sort(t)))/2

  ## using the median value of t as the initial value of a.
  a.old <- t[which(x == max(x))[1]]
  a.new <- a.old
  change <- 10 * epsilon

  ## n.iter is the number of iteration covered so far.
  n.iter <- 0

  while ((change > epsilon) & (n.iter < max.iter)) {
    a.old <- a.new
    n.iter <- n.iter + 1
    sigma <- solve.sigma(x, t, a.old)
    if (n.iter == 1) {
        sigma[is.na(sigma)] <- as.numeric(sigma[which(!is.na(sigma))])[1] / 10
    }
    a.new <- solve.a(x, t, a.old, sigma$sigma.1, sigma$sigma.2)
    change <- abs(a.old - a.new)
  }
  d <- x
  sigma$sigma.2 <- sqrt(sigma$sigma.2)
  sigma$sigma.1 <- sqrt(sigma$sigma.1)

  d[t < a.new] <- dnorm(t[t < a.new], mean = a.new, sd = sigma$sigma.1) * sigma$sigma.1
  d[t >= a.new] <- dnorm(t[t >= a.new], mean = a.new, sd = sigma$sigma.2) * sigma$sigma.2
  scale <- exp(sum(d[d > 1e-3]^2 * log(x[d > 1e-3] / d[d > 1e-3])) / sum(d[d > 1e-3]^2))
  return(c(a.new, sigma$sigma.1, sigma$sigma.2, scale))
}

#' Computes vector of cumulative sums on reversed input. Returns cumulative sum vector going from the sum of all elements to one.
#' @param x float - vector of numerical values
#' @return Returns a vector
#' @export
rev_cum_sum <- function(x) {
  x <- rev(x)
  return(rev(cumsum(x)))
}

#' TODO: Document
#' @export
compute_start_bound <- function(x, left_sigma_ratio_lim) {
  start_bound <- 1
  
  len_x <- length(x)
  idx <- which(x >= left_sigma_ratio_lim / (left_sigma_ratio_lim + 1) * x[len_x])
  if (length(idx) > 0) {
    start_bound <- max(1, min(idx))
  }
  return (start_bound)
}

#' TODO: Document
#' @export
compute_end_bound <- function(x, right_sigma_ratio_lim) {
  len_x <- length(x)
  end_bound <- len_x - 1

  idx <- which(x <= right_sigma_ratio_lim / (right_sigma_ratio_lim + 1) * x[len_x])
  if (length(idx) > 0) {
    end_bound <- min(len_x - 1, max(idx))
  }
  return (end_bound)
}

#' @param x Cumulative intensity values.
#' @param sigma.ratio.lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation.
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @return Returns a list with bounds with following items:
#' \itemize{
#'   \item start - start bound
#'   \item end - end bound
#'}
#' @export
compute_bounds <- function(x, sigma.ratio.lim) {
  start <- compute_start_bound(x, sigma.ratio.lim[1])
  end <- compute_end_bound(x, sigma.ratio.lim[2])
  return(list(start = start, end = end))
}

#' Compute difference between neighbouring elements of a vector and optionally apply a mask such that the maximum difference is no higher than 4-fold minimum difference.
#' @param x - float - a vector of numerical values.
#' @param apply_mask - boolean - whether to apply threshold mask to the output vector.
#' @return Returns vector of numeric differences between neighbouring values.
#' @export
compute_dx <- function(x, apply_mask=TRUE) {
  l <- length(x)
  diff_x <- diff(x)
  if (l == 2) {
      dx <- rep(diff_x, 2)
  } else {
    dx <- c(
      x[2] - x[1],
      diff(x, lag = 2) / 2,
      x[l] - x[l - 1]
    )
  }
  if (apply_mask) {
    diff_threshold <- min(diff_x) * 4
    dx <- pmin(dx, diff_threshold)
  }
  return (dx)
}

#' Find base.curve RTs that lay within RT range of the whole feature table and append the intensities to these RTs.
#' @param profile Profile table with shape number-of-features*4. The table contains following columns:
#' \itemize{
#'   \item mz - float - mass-to-charge ratio of feature
#'   \item rt - float - retention time of features
#'   \item intensity - float - intensity of features
#'   \item group_number - integer - group number assigned to each feature based on their rt similarity
#' }
#' @param base.curve Matrix that contains rts of feature in the same rt cluster.
#' @return dataframe with two columns
#' @export
compute_chromatographic_profile <- function(profile, base.curve) {
  rt_range <- range(profile[, "rt"])
  rt_profile <- base.curve[between(base.curve[, "base.curve"], min(rt_range), max(rt_range)), ]
  rt_profile[rt_profile[, "base.curve"] %in% profile[, "rt"], 2] <- profile[, "intensity"]
  colnames(rt_profile)[2] <- "intensity"

  return (rt_profile)
}

#' Estimate total signal strength (total area of the estimated normal curve).
#' @param y - float - a vector of intensities.
#' @param d - float - a vector of \emph{y} values in a gaussian curve.
#' @param scale - float - a vector of scaled intensity values.
#' @export
compute_scale <- function(y, d) {
  dy_ratio <- d^2 * log(y / d)
  dy_ratio[is.na(dy_ratio)] <- 0
  dy_ratio[is.infinite(dy_ratio)] <- 0

  scale <- exp(sum(dy_ratio) / sum(d^2))
  return (scale)
}

#' Estimate the parameters of Bi-Gaussian curve.
#' @param x Vector of RTs that lay in the same RT cluster.
#' @param y Intensities that belong to x.
#' @param power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param sigma.ratio.lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @return A vector with length 4. The items are as follows going from first to last:
#' \itemize{
#'   \item mean of gaussian curve
#'   \item standard deviation at the left side of the gaussian curve
#'   \item standard deviation at the right side of the gaussian curve
#'   \item estimated total signal strength (total area of the estimated normal curve)
#'}
#' @export
bigauss.esti <- function(x, y, power = 1, do.plot = FALSE, sigma.ratio.lim = c(0.3, 3)) {
  # even producing a dataframe with x and y as columns without actually using it causes the test to run forever
  sel <- which(y > 1e-10)
  if (length(sel) < 2) {
    return (c(median(x), 1, 1, 0))
  } else {
    x <- x[sel]
    y <- y[sel]

    y.0 <- y
    if (do.plot) {
        plot(x, y)
    }
    max.y.0 <- max(y.0, na.rm = TRUE)
    y <- (y / max.y.0)^power

    dx <- compute_dx(x)

    y.cum <- cumsum(y * dx)
    x.y.cum <- cumsum(y * x * dx)
    xsqr.y.cum <- cumsum(y * x^2 * dx)

    y.cum.rev <- rev_cum_sum(y * dx)
    x.y.cum.rev <- rev_cum_sum(x * y * dx)
    xsqr.y.cum.rev <- rev_cum_sum(y * x^2 * dx)

    bounds <- compute_bounds(y.cum, sigma.ratio.lim)
    end <- bounds$end
    start <- bounds$start

    if (end <= start) {
      m <- min(mean(x[start:end]), x[max(which(y.cum.rev > 0))])
    } else {
      m.candi <- x[start:end] + diff(x[start:(end + 1)]) / 2
      rec <- matrix(0, ncol = 3, nrow = end - start + 1)

      s1 <- sqrt((xsqr.y.cum[start:end] + m.candi^2 * y.cum[start:end] - 2 * m.candi * x.y.cum[start:end]) / y.cum[start:end])
      s2 <- sqrt((xsqr.y.cum.rev[start:end + 1] + m.candi^2 * y.cum.rev[start:end + 1] - 2 * m.candi * x.y.cum.rev[start:end + 1]) / y.cum.rev[start:end + 1])
      rec[, 1] <- s1
      rec[, 2] <- s2
      rec[, 3] <- y.cum[start:end] / y.cum.rev[start:end + 1]

      d <- log(rec[, 1] / rec[, 2]) - log(rec[, 3])
      if (min(d, na.rm = TRUE) * max(d, na.rm = TRUE) < 0) {
        sel <- c(which(d == max(d[d < 0]))[1], which(d == min(d[d >= 0])))
        m <- (sum(abs(d[sel]) * m.candi[sel])) / (sum(abs(d[sel])))
      } else {
        d <- abs(d)
        m <- m.candi[which(d == min(d, na.rm = TRUE))[1]]
      }
    }

    if (do.plot) {
        abline(v = m)
    }

    sel1 <- which(x < m)
    sel2 <- which(x >= m)
    s1 <- sqrt(sum((x[sel1] - m)^2 * y[sel1] * dx[sel1]) / sum(y[sel1] * dx[sel1]))
    s2 <- sqrt(sum((x[sel2] - m)^2 * y[sel2] * dx[sel2]) / sum(y[sel2] * dx[sel2]))

    s1 <- s1 * sqrt(power)
    s2 <- s2 * sqrt(power)

    d1 <- dnorm(x[sel1], sd = s1, mean = m)
    d2 <- dnorm(x[sel2], sd = s2, mean = m)
    d <- c(d1 * s1, d2 * s2) # notice this "density" doesnt integrate to 1. Rather it integrates to (s1+s2)/2
    y <- y.0

    scale <- compute_scale(y, d)

    if (do.plot) {
      lines(x[y > 0], d * scale, col = "red")
    }

    to.return <- c(m, s1, s2, scale)
    if (sum(is.na(to.return)) > 0) {
      m <- sum(x * y) / sum(y)
      s1 <- s2 <- sum(y * (x - m)^2) / sum(y)
      scale <- sum(y) / s1
      to.return <- c(m, s1, s2, scale)
    }
  }
  return(to.return)
}

#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param vlys A vector of sorted RT-valley values at which the kernel estimate was computed.
#' @param dx Difference between neighbouring RT values with step 2.
#' @param pks A vector of sorted RT-peak values at which the kernel estimate was computed.
#' @export
compute_initiation_params <- function(rt_profile, vlys, dx, pks) {
  m <- s1 <- s2 <- delta <- pks
  for (i in 1:length(m))
  {
    sel.1 <- which(rt_profile[, "base.curve"] >= max(vlys[vlys < m[i]]) & rt_profile[, "base.curve"] < m[i])
    s1[i] <- sqrt(sum((rt_profile[sel.1, "base.curve"] - m[i])^2 * rt_profile[sel.1, "intensity"] * dx[sel.1]) / sum(rt_profile[sel.1, "intensity"] * dx[sel.1]))

    sel.2 <- which(rt_profile[, "base.curve"] >= m[i] & rt_profile[, "base.curve"] < min(vlys[vlys > m[i]]))
    s2[i] <- sqrt(sum((rt_profile[sel.2, "base.curve"] - m[i])^2 * rt_profile[sel.2, "intensity"] * dx[sel.2]) / sum(rt_profile[sel.2, "intensity"] * dx[sel.2]))

    delta[i] <- (sum(rt_profile[sel.1, "intensity"] * dx[sel.1]) + sum(rt_profile[sel.2, "intensity"] * dx[sel.2])) / ((sum(dnorm(rt_profile[sel.1, "base.curve"], mean = m[i], sd = s1[i])) * s1[i] / 2) + (sum(dnorm(rt_profile[sel.2, "base.curve"], mean = m[i], sd = s2[i])) * s2[i] / 2))
  }
  return (list(s1 = s1,
    s2 = s2,
    delta = delta))
}

#' @param m A vector of sorted RT-peak values at which the kernel estimate was computed.
#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param delta Parameter computed by the initiation step.
#' @param s1 Parameter computed by the initiation step.
#' @param s2 Parameter computed by the initiation step.
#' @export
compute_e_step <- function(m, rt_profile, delta, s1, s2) {
  fit <- matrix(0, ncol = length(m), nrow = length(rt_profile[, "base.curve"])) # this is the matrix of fitted values
  cuts <- c(-Inf, m, Inf)
  for (j in 2:length(cuts))
  {
    sel <- which(rt_profile[, "base.curve"] >= cuts[j - 1] & rt_profile[, "base.curve"] < cuts[j])
    use.s1 <- which(1:length(m) >= (j - 1))
    s.to.use <- s2
    s.to.use[use.s1] <- s1[use.s1]
    for (i in 1:ncol(fit))
    {
      fit[sel, i] <- dnorm(rt_profile[sel, "base.curve"], mean = m[i], sd = s.to.use[i]) * s.to.use[i] * delta[i]
    }
  }
  fit[is.na(fit)] <- 0
  sum.fit <- apply(fit, 1, sum)
  return(list(fit = fit, sum.fit = sum.fit))
}

#' @param rt_profile Dataframe that stores RTs and intensities of features.
#' @param power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param sigma.ratio.lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @param bw Bandwidth vector to use in the kernel smoother.
#' @param eliminate When a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param max.iter Maximum number of iterations when executing the E step.
#' @param estim.method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param BIC.factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#'  models with more peaks are penalized more.
#' @importFrom dplyr filter arrange
#' @export
bigauss.mix <- function(rt_profile, power = 1, do.plot = FALSE, sigma.ratio.lim = c(0.1, 10), bw = c(15, 30, 60), eliminate = .05, max.iter = 25, estim.method, BIC.factor = 2) {
  all.bw <- sort(bw)
  results <- new("list")
  smoother.pk.rec <- smoother.vly.rec <- new("list")
  bic.rec <- all.bw

  if (do.plot) {
    par(mfrow = c(ceiling(length(all.bw) / 2), 2))
    par(mar = c(1, 1, 1, 1))
  }

  last.num.pks <- Inf

  rt_profile_unfiltered <- rt_profile
  rt_profile <- data.frame(rt_profile) |>
    filter(intensity > 1e-5) |>
    arrange(base.curve)

  for (bw.n in length(all.bw):1)
  {
    bw <- all.bw[bw.n]
    this.smooth <- ksmooth(rt_profile_unfiltered[, "base.curve"], rt_profile_unfiltered[, "intensity"], kernel = "normal", bandwidth = bw)
    turns <- find.turn.point(this.smooth$y)
    pks <- this.smooth$x[turns$pks]
    vlys <- c(-Inf, this.smooth$x[turns$vlys], Inf)

    smoother.pk.rec[[bw.n]] <- pks
    smoother.vly.rec[[bw.n]] <- vlys
    if (length(pks) != last.num.pks) {
      last.num.pks <- length(pks)
      l <- length(rt_profile[, "base.curve"])
      dx <- compute_dx(rt_profile[, "base.curve"], apply_mask = FALSE)

      # initiation
      initiation_params <- compute_initiation_params(rt_profile, vlys, dx, pks)
      s1 <- initiation_params$s1
      s2 <- initiation_params$s2
      delta <- initiation_params$delta

      delta[is.na(delta)] <- 1e-10
      s1[is.na(s1)] <- 1e-10
      s2[is.na(s2)] <- 1e-10

      this.change <- Inf
      counter <- 0

      m <- pks
      while (this.change > 0.1 & counter <= max.iter) {
        counter <- counter + 1
        old.m <- m

        # E step
        fits <- compute_e_step(m, rt_profile, delta, s1, s2)
        fit <- fits$fit
        sum.fit <- fits$sum.fit

        # Elimination step
        fit <- fit / sum.fit
        fit2 <- fit * rt_profile[, "intensity"]
        perc.explained <- apply(fit2, 2, sum) / sum(rt_profile[, "intensity"])
        max.erase <- max(1, round(length(perc.explained) / 5))

        to.erase <- which(perc.explained <= min(eliminate, perc.explained[order(perc.explained, na.last = FALSE)[max.erase]]))


        if (length(to.erase) > 0) {
          m <- m[-to.erase]
          s1 <- s1[-to.erase]
          s2 <- s2[-to.erase]
          delta <- delta[-to.erase]
          fit <- fit[, -to.erase]
          if (is.null(ncol(fit))) {
            fit <- matrix(fit, ncol = 1)
          }
          sum.fit <- apply(fit, 1, sum)
          fit <- fit / sum.fit
          old.m <- old.m[-to.erase]
        }

        # M step
        for (i in 1:length(m))
        {
          this.y <- rt_profile[, "intensity"] * fit[, i]
          if (estim.method == "moment") {
            this.fit <- bigauss.esti(rt_profile[, "base.curve"], this.y, power = power, do.plot = FALSE, sigma.ratio.lim = sigma.ratio.lim)
          } else {
            this.fit <- bigauss.esti.EM(rt_profile[, "base.curve"], this.y, power = power, do.plot = FALSE, sigma.ratio.lim = sigma.ratio.lim)
          }
          m[i] <- this.fit[1]
          s1[i] <- this.fit[2]
          s2[i] <- this.fit[3]
          delta[i] <- this.fit[4]
        }
        delta[is.na(delta)] <- 0

        # amount of change
        this.change <- sum((old.m - m)^2)
      }
      cuts <- c(-Inf, m, Inf)
      fit <- fit * 0
      for (j in 2:length(cuts))
      {
        sel <- which(rt_profile[, "base.curve"] >= cuts[j - 1] & rt_profile[, "base.curve"] < cuts[j])
        use.s1 <- which(1:length(m) >= (j - 1))
        s.to.use <- s2
        s.to.use[use.s1] <- s1[use.s1]
        for (i in 1:ncol(fit))
        {
          if (s.to.use[i] != 0) {
            fit[sel, i] <- dnorm(rt_profile[sel, "base.curve"], mean = m[i], sd = s.to.use[i]) * s.to.use[i] * delta[i]
          }
        }
      }

      if (do.plot) {
        plot_rt_profile(rt_profile, bw, fit, m)
      }
      area <- delta * (s1 + s2) / 2
      rss <- sum((rt_profile[, "intensity"] - apply(fit, 1, sum))^2)
      l <- length(rt_profile[, "base.curve"])
      bic <- l * log(rss / l) + 4 * length(m) * log(l) * BIC.factor
      results[[bw.n]] <- cbind(m, s1, s2, delta, area)
      bic.rec[bw.n] <- bic
    } else {
      results[[bw.n]] <- NA
      bic.rec[bw.n] <- Inf
      results[[bw.n]] <- results[[bw.n + 1]]
    }
  }
  sel <- which(bic.rec == min(bic.rec, na.rm = TRUE))
  if (length(sel) > 1) {
    sel <- sel[which(all.bw[sel] == max(all.bw[sel]))]
  }
  rec <- new("list")
  rec$param <- results[[sel]]
  rec$smoother.pks <- smoother.pk.rec
  rec$smoother.vlys <- smoother.vly.rec
  rec$all.param <- results
  rec$bic <- bic.rec
  return(rec)
}

#' Reevaluate parameters of rtomatographic gaussian curves.
#' @param that.curve Dataframe that stores RTs and intensities of features.
#' @param pks A vector of sorted RT-peak values at which the kernel estimate was computed.
#' @param vlys A vector of sorted RT-valley values at which the kernel estimate was computed.
#' @param ignore In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a
#' proportion of intensities less than this value, the component will be ignored.
#' @param max.iter Maximum number of iterations when reevaluating gaussian curves.
#' @param aver_diff Average retention time difference across RTs of all features.
#' @importFrom dplyr between
#' @export
normix <- function(that.curve, pks, vlys, ignore = 0.1, max.iter = 50, aver_diff) {
  x <- that.curve[, 1]
  y <- that.curve[, 2]
  rt_int_list <- list(rt = x, intensities = y)

  if (length(pks) == 1) {
    mu_sc_std <- compute_mu_sc_std(rt_int_list, aver_diff)
    miu <- mu_sc_std$label
    sc <- mu_sc_std$intensity
    sigma <- mu_sc_std$sigma
  } else {
    pks <- sort(pks)
    vlys <- sort(vlys)
    l <- length(pks)
    miu <- sigma <- sc <- pks
    w <- matrix(0, nrow = l, ncol = length(x))


    for (m in 1:l)
    {
      # this pattern occurs multiple times in other scripts
      this.low <- max(vlys[vlys <= pks[m]])
      this.high <- min(vlys[vlys >= pks[m]])

      indices <- dplyr::between(x, this.low, this.high)
      this.x <- x[indices]
      this.y <- y[indices]

      if (length(this.x) == 0 | length(this.y) == 0) {
        miu[m] <- NaN
        sigma[m] <- NaN
        sc[m] <- 1
      } else {
        rt_int_list_this <- list(rt = this.x, intensities = this.y)
        mu_sc_std <- compute_mu_sc_std(rt_int_list_this, aver_diff)
        miu[m] <- mu_sc_std$label
        sc[m] <- mu_sc_std$intensity
        sigma[m] <- mu_sc_std$sigma
      }
    }

    to.erase <- which(is.na(miu) | is.na(sigma) | sigma == 0 | is.na(sc))
    if (length(to.erase) > 0) {
      l <- l - length(to.erase)
      miu <- miu[-to.erase]
      sigma <- sigma[-to.erase]
      sc <- sc[-to.erase]
      w <- w[-to.erase, ]
    }

    direc <- 1
    diff <- 1000
    iter <- 0

    while (diff > 0.05 & iter < max.iter) {
      iter <- iter + 1
      if (l == 1) {
        mu_sc_std <- compute_mu_sc_std(rt_int_list, aver_diff)
        miu <- mu_sc_std$label
        sc <- mu_sc_std$intensity
        sigma <- mu_sc_std$sigma
        break
      }
      miu.0 <- miu
      sigma.0 <- sigma
      sc.0 <- sc

      all.w <- y * 0
      for (m in 1:l)
      {
        all.w <- all.w + dnorm(x, mean = miu[m], sd = sigma[m]) * sc[m]
      }

      # when l is zero the iteration goes from 1 to 0 znd results in "index out of bound" error
      for (m in 1:l)
      {
        w[m, ] <- dnorm(x, mean = miu[m], sd = sigma[m]) * sc[m] / all.w
      }

      if (sum(is.na(w)) > 0) {
        break
      }

      for (m in 1:l)
      {
        this.y <- y * w[m, ]
        rt_int_list_this <- list(rt = x, intensities = this.y)
        mu_sc_std <- compute_mu_sc_std(rt_int_list_this, aver_diff)
        miu[m] <- mu_sc_std$label
        sc[m] <- mu_sc_std$intensity
        sigma[m] <- mu_sc_std$sigma

        if (sigma[m] == 0) {
          sc[m] <- NA
        }
      }
      diff <- sum((miu.0 - miu)^2)

      www <- w
      for (m in 1:l)
      {
        www[m, ] <- www[m, ] * y
      }
      www <- apply(www, 1, sum)
      www[which(is.na(sc))] <- 0
      www <- www / sum(www)
      max.erase <- max(1, round(l / 5))

      to.erase <- which(www <= min(ignore, www[order(www, na.last = FALSE)[max.erase]]))

      if (length(to.erase) > 0) {
        l <- l - length(to.erase)
        miu <- miu[-to.erase]
        sigma <- sigma[-to.erase]
        sc <- sc[-to.erase]
        w <- w[-to.erase, ]
        diff <- 1000
      }
    }
  }
  l <- length(miu)
  if (l == 1) {
    rec <- matrix(c(miu, sigma, sc), nrow = 1)
  } else {
    rec <- cbind(miu, sigma, sc)
  }
  colnames(rec) <- c("miu", "sigma", "scale")
  return(rec)
}

#' @param x Vector of RTs that lay in the same RT cluster.
#' @param y Intensities that belong to x.
#' @param power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param bw Bandwidth vector to use in the kernel smoother.
#' @param eliminate When a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param max.iter Maximum number of iterations when executing the E step.
#' @param BIC.factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#' @param aver_diff Average retention time difference across RTs of all features.
#' @export
normix.bic <- function(x, y, power = 2, do.plot = FALSE, bw = c(15, 30, 60), eliminate = .05, max.iter = 50, BIC.factor = 2, aver_diff) {
  all.bw <- bw[order(bw)]
  sel <- y > 1e-5
  x <- x[sel]
  y <- y[sel]
  sel <- order(x)
  y <- y[sel]
  x <- x[sel]
  results <- new("list")
  smoother.pk.rec <- smoother.vly.rec <- new("list")
  bic.rec <- all.bw

  if (do.plot) {
    par(mfrow = c(ceiling(length(all.bw) / 2), 2))
    par(mar = c(1, 1, 1, 1))
  }

  last.num.pks <- Inf

  for (bw.n in length(all.bw):1)
  {
    bw <- all.bw[bw.n]
    this.smooth <- ksmooth(x, y, kernel = "normal", bandwidth = bw)
    turns <- find.turn.point(this.smooth$y)
    pks <- this.smooth$x[turns$pks]
    vlys <- c(-Inf, this.smooth$x[turns$vlys], Inf)

    smoother.pk.rec[[bw.n]] <- pks
    smoother.vly.rec[[bw.n]] <- vlys
    if (length(pks) != last.num.pks) {
      last.num.pks <- length(pks)
      aaa <- normix(cbind(x, y), pks = pks, vlys = vlys, ignore = eliminate, max.iter = max.iter, aver_diff = aver_diff)

      total.fit <- x * 0
      for (i in 1:nrow(aaa))
      {
        total.fit <- total.fit + dnorm(x, mean = aaa[i, 1], sd = aaa[i, 2]) * aaa[i, 3]
      }

      if (do.plot) {
        plot_normix_bic(x, y, bw, aaa)
      }

      rss <- sum((y - total.fit)^2)
      l <- length(x)
      bic <- l * log(rss / l) + 3 * nrow(aaa) * log(l) * BIC.factor
      results[[bw.n]] <- aaa
      bic.rec[bw.n] <- bic
    } else {
      bic.rec[bw.n] <- Inf
      results[[bw.n]] <- results[[bw.n + 1]]
    }
  }
  sel <- which(bic.rec == min(bic.rec))
  if (length(sel) > 1) {
    sel <- sel[which(all.bw[sel] == max(all.bw[sel]))]
  }
  rec <- new("list")
  rec$param <- results[[sel]]
  rec$smoother.pks <- smoother.pk.rec
  rec$smoother.vlys <- smoother.vly.rec
  rec$all.param <- results
  rec$bic <- bic.rec
  return(rec)
}

#' Generate feature table from noise-removed LC/MS profile.
#'
#' @description
#' Each LC/MS profile is first processed by the function proc.cdf() to remove noise and reduce data size. A matrix containing m/z
#' value, retention time, intensity, and group number is output from proc.cdf(). This matrix is then fed to the function
#' prof.to.features() to generate a feature list. Every detected feature is summarized into a single row in the output matrix from this function.
#'
#' @param profile The matrix output from proc.cdf(). It contains columns of m/z value, retention time, intensity and group number.
#' @param bandwidth A value between zero and one. Multiplying this value to the length of the signal along the time axis helps
#'  determine the bandwidth in the kernel smoother used for peak identification.
#' @param min.bw The minimum bandwidth to use in the kernel smoother.
#' @param max.bw The maximum bandwidth to use in the kernel smoother.
#' @param sd.cut A vector of two. Features with standard deviation outside the range defined by the two numbers are eliminated.
#' @param sigma.ratio.lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @param shape.model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian".
#'  When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param estim.method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param do.plot Whether to generate diagnostic plots.
#' @param power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param component.eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a
#'  proportion of intensities less than this value, the component will be ignored.
#' @param BIC.factor the factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#'  models with more peaks are penalized more.
#' @return A matrix is returned. The columns are: m/z value, retention time, spread (standard deviation of the estimated normal
#'  curve), and estimated total signal strength (total area of the estimated normal curve).
#' @export
#' @examples
#' prof.to.features(extracted_features, sd.cut = sd_cut, sigma.ratio.lim = sigma_ratio_lim, do.plot = FALSE)
prof.to.features <- function(profile,
                             bandwidth = 0.5,
                             min.bw = NA,
                             max.bw = NA,
                             sd.cut = c(0.01, 500),
                             sigma.ratio.lim = c(0.01, 100),
                             shape.model = "bi-Gaussian",
                             estim.method = "moment",
                             do.plot = TRUE,
                             power = 1,
                             component.eliminate = 0.01,
                             BIC.factor = 2) {
  validate_inputs(shape.model, estim.method)

  profile <- preprocess_profile(profile)

  bws <- preprocess_bandwidth(min.bw, max.bw, profile)
  min.bw <- bws[["min.bw"]]
  max.bw <- bws[["max.bw"]]

  # base.curve <- compute_base_curve(profile[, "rt"])
  base.curve <- sort(unique(profile[, "rt"]))
  base.curve <- cbind(base.curve, base.curve * 0)
  all_rts <- compute_delta_rt(base.curve[, 1])
  aver_diff <- mean(diff(base.curve))

  keys <- c("mz", "pos", "sd1", "sd2", "area")
  peak_parameters <- matrix(0, nrow = 0, ncol = length(keys), dimnames = list(NULL, keys))

  feature_groups <- split(profile, profile$group_number)
  for (i in seq_along(feature_groups))
  {
    feature_group <- feature_groups[[i]]
    feature_group <- feature_group[order(feature_group[, "rt"]), ]

    num_features <- nrow(feature_group)
    if (dplyr::between(num_features, 2, 10)) {
      eic_area <- interpol.area(feature_group[, "rt"], feature_group[, "intensity"], base.curve[, "base.curve"], all_rts)
      rt_peak_shape <- c(median(feature_group[, "mz"]), median(feature_group[, "rt"]), sd(feature_group[, "rt"]), sd(feature_group[, "rt"]), eic_area)
      peak_parameters <- rbind(peak_parameters, rt_peak_shape)
    }
    if (num_features < 2) {
      time_weights <- all_rts[which(base.curve[, "base.curve"] %in% feature_group[2])]
      rt_peak_shape <- c(feature_group[1], feature_group[2], NA, NA, feature_group[3] * time_weights)
      peak_parameters <- rbind(peak_parameters, rt_peak_shape)
    }
    if (num_features > 10) {
      rt_range <- range(feature_group[, "rt"])
      bw <- min(max(bandwidth * (max(rt_range) - min(rt_range)), min.bw), max.bw)
      bw <- seq(bw, 2 * bw, length.out = 3)
      if (bw[1] > 1.5 * min.bw) {
        bw <- c(max(min.bw, bw[1] / 2), bw)
      }

      rt_profile <- compute_chromatographic_profile(feature_group, base.curve)
      if (shape.model == "Gaussian") {
        rt_peak_shape <- compute_gaussian_peak_shape(rt_profile, power, bw, component.eliminate, BIC.factor, aver_diff)
      } else {
        rt_peak_shape <- bigauss.mix(rt_profile, sigma.ratio.lim = sigma.ratio.lim, bw = bw, power = power, estim.method = estim.method, eliminate = component.eliminate, BIC.factor = BIC.factor)$param[, c(1, 2, 3, 5)]
      }

      if (is.null(nrow(rt_peak_shape))) {
        peak_parameters <- rbind(peak_parameters, c(median(feature_group[, "mz"]), rt_peak_shape))
      } else {
        for (m in 1:nrow(rt_peak_shape))
        {
          rt_diff <- abs(feature_group[, "rt"] - rt_peak_shape[m, 1])
          peak_parameters <- rbind(peak_parameters, c(mean(feature_group[which(rt_diff == min(rt_diff)), 1]), rt_peak_shape[m, ]))
        }
      }
    }
  }
  peak_parameters <- peak_parameters[order(peak_parameters[, "mz"], peak_parameters[, "pos"]), ]
  peak_parameters <- peak_parameters[which(apply(peak_parameters[, c("sd1", "sd2")], 1, min) > sd.cut[1] & apply(peak_parameters[, c("sd1", "sd2")], 1, max) < sd.cut[2]), ]
  rownames(peak_parameters) <- NULL

  if (do.plot) {
    plot_peak_summary(feature_groups, peak_parameters)
  }

  return(tibble::as_tibble(peak_parameters))
}
