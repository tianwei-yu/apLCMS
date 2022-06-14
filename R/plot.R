plot_raw_profile_histogram <- function(
  raw.prof,
  min.pres,
  baseline.correct,
  baseline.correct.noise.percentile,
  tol,
  new.prof
) {
  h.1 <- log10(raw.prof$height.rec[raw.prof$height.rec[, 2] <= max(2, raw.prof$min.count.run * min.pres / 2), 3])
  h.2 <- log10(raw.prof$height.rec[raw.prof$height.rec[, 2] >= raw.prof$min.count.run * min.pres, 3])

  if (is.na(baseline.correct)) {
    baseline.correct <- 10^quantile(h.1, baseline.correct.noise.percentile)
    message(c("maximal height cut is automatically set at the", baseline.correct.noise.percentile, "percentile of noise group heights: ", baseline.correct))
  } else {
    message(c("maximal height cut is provided by user: ", baseline.correct))
  }
  par(mfrow = c(2, 2))

  plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "tolerance level loaded", axes = FALSE)
  text(x = 0, y = 0, tol, cex = 1.2)

  if (length(h.1) > 50) {
    plot(density(h.1), xlab = "maximum height of group (log scale)", xlim = range(c(h.1, h.2)), main = "Black - noise groups \n Blue - selected groups")
  } else {
    plot(NA, NA, xlab = "maximum height of group (log scale)", xlim = range(c(h.1, h.2)), ylim = c(0, 1), main = "Black - noise groups \n Blue - selected groups")
    if (length(h.1) > 0) abline(v = h.1)
  }

  abline(v = log10(baseline.correct), col = "red")
  lines(density(log10(new.prof$height.rec)), col = "blue")
  hist(new.prof$time.range.rec, xlab = "Range of retention time in the same group", ylab = "Density", freq = FALSE, nclass = 100, main = "Group retention time range distribution")
  hist(new.prof$mz.pres.rec, xlab = "% signal present in the same group", ylab = "Density", freq = FALSE, nclass = 20, main = "Group % present signal distribution")
}