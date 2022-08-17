draw_plot <- function(x = c(-1, 1), y = c(-1, 1),
                      xlab = "", ylab = "",
                      main = "", axes = FALSE,
                      type = "n", label = NA, cex = 1.2) {
  plot(x, y, type = type, xlab = xlab, ylab = ylab, main = main, axes = axes)
  if (!is.na(label)) {
    text(x = 0, y = 0, label, cex = cex)
  }
}

tolerance_plot <- function(x, y, exp_y, selected, main) {
  plot(x, y, xlab = "Delta", ylab = "Density", main = main, cex = .25)
  lines(x, exp_y, col = "red")
  abline(v = x[selected], col = "blue")
}

#' @export
draw_chr_normal_peaks <- function(x, truth) {
  true.y1 <- dnorm(x[x < truth[1]], mean = truth[1], sd = truth[2]) * truth[2] * truth[4]
  true.y2 <- dnorm(x[x >= truth[1]], mean = truth[1], sd = truth[3]) * truth[3] * truth[4]
  lines(x, c(true.y1, true.y2), col = "green")
}

#' @export
plot_raw_profile_histogram <- function(raw.prof,
                                       min.pres,
                                       baseline.correct,
                                       baseline.correct.noise.percentile,
                                       tol,
                                       new.prof) {
  h.1 <- log10(raw.prof$height.rec[raw.prof$height.rec[, 2] <= max(2, raw.prof$min.count.run * min.pres / 2), 3])
  h.2 <- log10(raw.prof$height.rec[raw.prof$height.rec[, 2] >= raw.prof$min.count.run * min.pres, 3])
  
  if (is.na(baseline.correct)) {
    baseline.correct <- 10 ^ quantile(h.1, baseline.correct.noise.percentile)
    message(c("maximal height cut is automatically set at the",
              baseline.correct.noise.percentile,
              "percentile of noise group heights: ",
              baseline.correct
             )
           )
  } else {
    message(c("maximal height cut is provided by user: ", baseline.correct))
  }
  par(mfrow = c(2, 2))
  
  draw_plot(main = "tolerance level loaded", label = tol)
  
  if (length(h.1) > 50) {
    plot(density(h.1),
         xlab = "maximum height of group (log scale)",
         xlim = range(c(h.1, h.2)),
         main = "Black - noise groups \n Blue - selected groups")
  } else {
    plot(NA, NA, xlab = "maximum height of group (log scale)",
      xlim = range(c(h.1, h.2)), ylim = c(0, 1),
      main = "Black - noise groups \n Blue - selected groups"
        )
    if (length(h.1) > 0)
      abline(v = h.1)
  }
  
  abline(v = log10(baseline.correct), col = "red")
  lines(density(log10(new.prof$height.rec)), col = "blue")
  hist(
    new.prof$time.range.rec,
    xlab = "Range of retention time in the same group",
    ylab = "Density",
    freq = FALSE,
    nclass = 100,
    main = "Group retention time range distribution"
  )
  hist(
    new.prof$mz.pres.rec,
    xlab = "% signal present in the same group",
    ylab = "Density",
    freq = FALSE,
    nclass = 20,
    main = "Group % present signal distribution"
  )
}

#' @export
plot_peak_summary <- function(feature_groups, processed_features) {
  mz_sd <- compute_mz_sd(feature_groups)

  par(mfrow = c(2, 2))
  plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
  text(x = 0, y = 0, "Estimate peak \n area/location", cex = 1.5)
  hist(mz_sd, xlab = "m/z SD", ylab = "Frequency", main = "m/z SD distribution")
  hist(c(processed_features[, "sd1"], processed_features[, "sd2"]), xlab = "Retention time SD", ylab = "Frequency", main = "Retention time SD distribution")
  hist(log10(processed_features[, "area"]), xlab = "peak strength (log scale)", ylab = "Frequency", main = "Peak strength distribution")
}

#' @export
plot_chr_profile <- function(chr_profile, bw, fit, m) {
  plot(chr_profile[, "base_curve"], chr_profile[, "intensity"], cex = .1, main = paste("bw=", bw))
  sum.fit <- apply(fit, 1, sum)
  lines(chr_profile[, "base_curve"], sum.fit)
  abline(v = m)
  cols <- c("red", "green", "blue", "cyan", "brown", "black", rep("grey", 100))
  for (i in 1:length(m))
  {
    lines(chr_profile[, "base_curve"], fit[, i], col = cols[i])
  }
}

#' @export
plot_normix_bic <- function(x, y, bw, aaa) {
  plot(x, y, cex = .1, main = paste("bw=", bw))
  abline(v = aaa[, 1])
  cols <- c("red", "green", "blue", "cyan", "brown", "black", rep("grey", 100))
  for (i in 1:nrow(aaa))
  {
    lines(x, dnorm(x, mean = aaa[i, 1], sd = aaa[i, 2]) * aaa[i, 3], col = cols[i])
  }
}
