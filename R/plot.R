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
