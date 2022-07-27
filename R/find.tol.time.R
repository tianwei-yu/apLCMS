#' An internal function that find elution time tolerance level. 
#' 
#' This function finds the time tolerance level. Also, it returns the grouping information given the time tolerance.
#' 
#' @param mz mz value of all peaks in all profiles in the study.
#' @param chr retention time of all peaks in all profiles in the study.
#' @param lab label of all peaks in all profiles in the study.
#' @param number_of_samples The number of spectra in this analysis.
#' @param mz_tol_relative m/z tolerance level for the grouping of signals into peaks. This value is expressed as the percentage of the m/z value. 
#'  This value, multiplied by the m/z value, becomes the cutoff level.
#' @param rt_tol_relative the elution time tolerance. If NA, the function finds the tolerance level first. If a numerical value is given, 
#'  the function directly goes to the second step - grouping peaks based on the tolerance.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins the maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#' @param mz_tol_absolute As the m/z tolerance in alignment is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. 
#'  This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.
#' @param max.num.segments the maximum number of segments.
#' @param do.plot Indicates whether plot should be drawn.
#' @return A matrix with six columns. Every row corrsponds to a peak in one of the spectrum. The columns are: m/z, elution time, spread, signal strength, 
#'  spectrum label, and peak group label. The rows are ordered by the median m/z of each peak group, and with each peak group the rows are ordered 
#'  by the elution time.
#' @examples
#' find.tol.time(mz, chr, lab, number_of_samples = number_of_samples, mz_tol_relative = mz_tol_relative, mz_tol_absolute = mz_tol_absolute, do.plot = FALSE)
find.tol.time <- function(mz,
                          chr,
                          lab,
                          number_of_samples,
                          mz_tol_relative = 2e-5,
                          rt_tol_relative = NA,
                          aver.bin.size = 200,
                          min.bins = 50,
                          max.bins = 100,
                          mz_tol_absolute = 0.01,
                          max.num.segments = 10000,
                          do.plot = TRUE) {
    o <- order(mz)
    mz <- mz[o]
    chr <- chr[o]
    lab <- lab[o]
    rm(o)
    
    l <- length(mz)
    
    breaks <-
        c(0, which((mz[2:l] - mz[1:(l - 1)]) > min(mz_tol_absolute, mz_tol_relative * ((
            mz[2:l] + mz[1:(l - 1)]
        ) / 2))), l)
    
    for (i in 2:length(breaks)) {
        this.o <- order(chr[(breaks[i - 1] + 1):breaks[i]])
        this.o <- this.o + breaks[i - 1]
        mz[(breaks[i - 1] + 1):breaks[i]] <- mz[this.o]
        chr[(breaks[i - 1] + 1):breaks[i]] <- chr[this.o]
        lab[(breaks[i - 1] + 1):breaks[i]] <- lab[this.o]
    }
    
    breaks <- breaks[c(-1, -length(breaks))]
    if (is.na(rt_tol_relative)) {
        da <- 0
        if (length(breaks) > max.num.segments) {
            s <- floor(seq(2, length(breaks), length.out = max.num.segments))
        } else {
            s <- 2:length(breaks)
        }
        
        for (i in s) {
            this.sel <- (breaks[i - 1] + 1):breaks[i]
            
            if (length(this.sel) <= 3 * number_of_samples) {
                this.d <- as.vector(dist(chr[this.sel]))
                if (length(this.d) > 100)
                    this.d <- sample(this.d, 100)
                da <- c(da, this.d)
            }
        }
        
        da <- da[!is.na(da)]
        uppermost <- max(da)
        n = min(max.bins, max(min.bins, round(length(da) / aver.bin.size)))
        des <- density(da, kernel = "gaussian", n = n,
                       bw = uppermost / n * 2, from = 0)
        y <- des$y[des$x > 0]
        x <- des$x[des$x > 0]
        
        this.l <- lm(y[x > uppermost / 4] ~ x[x > uppermost / 4])
        
        exp.y <- this.l$coef[1] + this.l$coef[2] * x
        
        y2 <- y[1:(length(y) - 1)]
        y3 <- y[2:(length(y))]
        y2[which(y2 < y3)] <- y3[which(y2 < y3)]
        y[1:(length(y) - 1)] <- y2
        
        yy <- cumsum(y > 1.5 * exp.y)
        yi <- seq_along(yy)
        sel <- min(which(yy < yi)) - 1
        rt_tol_relative <- x[sel]
        
        if (do.plot) {
          tolerance_plot(x, y, exp.y, sel, main = "find retention time tolerance")
        }
    }
    
    da <- chr[2:l] - chr[1:(l - 1)]
    breaks.2 <- which(da > rt_tol_relative)
    all.breaks <- c(0, unique(c(breaks, breaks.2)), l)
    all.breaks <- all.breaks[order(all.breaks)]
    
    grps <- seq_along(mz)
    for (i in 2:length(all.breaks)) {
        grps[(all.breaks[i - 1] + 1):all.breaks[i]] <- i
    }
    
    list(
        mz = mz,
        chr = chr,
        lab = lab,
        grps = grps,
        chr.tol = rt_tol_relative,
        mz.tol = mz_tol_relative
    )
}
