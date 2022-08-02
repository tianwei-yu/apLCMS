
validate_inputs <- function(shape.model, estim.method) {
    if (!shape.model %in% c("Gaussian", "bi-Gaussian")) {
        stop("shape.model argument must be 'Gaussian' or 'bi-Gaussian'")
    }
    if (!estim.method %in% c("moment", "EM")) {
        stop("estim.method argument must be 'moment' or 'EM'")
    }
}

preprocess_bandwidth <- function(min.bw, max.bw, feature_table) {
    if (is.na(min.bw)) {
        min.bw <- diff(range(feature_table[, 2], na.rm = TRUE)) / 60
    }
    if (is.na(max.bw)) {
        max.bw <- diff(range(feature_table[, 2], na.rm = TRUE)) / 15
    }
    if (min.bw >= max.bw) {
        min.bw <- max.bw / 4
    }

    return (list("min.bw" = min.bw, "max.bw" = max.bw))
}

preprocess_feature_table <- function(feature_table) {
    keys <- c("mz", "rt", "intensity", "group_number")
    colnames(feature_table) <- keys

    return (data.frame(feature_table))
}

compute_gaussian_peak_shape <- function(this.curve, power, bw, component.eliminate, BIC.factor) {

    ## this function computes parameters of chromatographic peak shape if peaks are considered to be gaussian

    chr_peak_shape <- normix.bic(this.curve[, "base_curve"], this.curve[, 2], power = power, bw = bw, eliminate = component.eliminate, BIC.factor = BIC.factor)$param
    if (nrow(chr_peak_shape) == 1) {
        chr_peak_shape <- c(chr_peak_shape[1, 1:2], chr_peak_shape[1, 2], chr_peak_shape[1, 3])
    } else {
        chr_peak_shape <- cbind(chr_peak_shape[, 1:2], chr_peak_shape[, 2], chr_peak_shape[, 3])
    }
    return (chr_peak_shape)
}

solve.a <- function(x, t, a, sigma.1, sigma.2) {
    ## thif function solves the value of a using the x, t, a from the
    ## previous step, and sigma.1, and sigma.2

    w <- x * (as.numeric(t < a) / sigma.1 + as.numeric(t >= a) / sigma.2)
    return(sum(t * w) / sum(w))
}

#' @export
prep.uv <- function(x, t, a) {

    ## this function prepares the parameters required for latter
    ## compuation. u, v, and sum of x.

    temp <- (t - a)^2 * x
    u <- sum(temp * as.numeric(t < a))
    v <- sum(temp * as.numeric(t >= a))
    return(list(
        u = u,
        v = v,
        x.sum = sum(x)
    ))
}

#' @export
solve.sigma <- function(x, t, a) {
    ## this function takes the value intensity level x, retention time t
    ## and assumed breaking point a, calculates the square estimated of
    ## sigma.1 and sigma.2.


    tt <- prep.uv(x, t, a)
    sigma.1 <- tt$u / tt$x.sum * ((tt$v / tt$u)^(1 / 3) + 1)
    sigma.2 <- tt$v / tt$x.sum * ((tt$u / tt$v)^(1 / 3) + 1)
    return(list(
        sigma.1 = sigma.1,
        sigma.2 = sigma.2
    ))
}

#' @export
bigauss.esti.EM <- function(t, x, max.iter = 50, epsilon = 0.005, power = 1, do.plot = FALSE, truth = NA, sigma.ratio.lim = c(0.3, 1)) {
    ## function takes into x and t, and then computes the value of
    ## sigma.1, sigma.2 and a using iterative method. the returned
    ## values include estimated sigmas, a and a boolean variable on
    ## whether the termination criteria is satified upon the end of the
    ## program.

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
        ## print(c(n.iter,change))
        a.old <- a.new
        n.iter <- n.iter + 1
        sigma <- solve.sigma(x, t, a.old)
        if (n.iter == 1) sigma[is.na(sigma)] <- as.numeric(sigma[which(!is.na(sigma))])[1] / 10
        a.new <- solve.a(x, t, a.old, sigma$sigma.1, sigma$sigma.2)
        change <- abs(a.old - a.new)
    }
    #  return(list(a=a.new,
    #              sigma.1=sigma$sigma.1,
    #              sigma.2=sigma$sigma.2,
    #              iter.end=(max.iter>n.iter)
    #              ))
    d <- x
    sigma$sigma.2 <- sqrt(sigma$sigma.2)
    sigma$sigma.1 <- sqrt(sigma$sigma.1)

    d[t < a.new] <- dnorm(t[t < a.new], mean = a.new, sd = sigma$sigma.1) * sigma$sigma.1
    d[t >= a.new] <- dnorm(t[t >= a.new], mean = a.new, sd = sigma$sigma.2) * sigma$sigma.2
    scale <- exp(sum(d[d > 1e-3]^2 * log(x[d > 1e-3] / d[d > 1e-3])) / sum(d[d > 1e-3]^2))
    return(c(a.new, sigma$sigma.1, sigma$sigma.2, scale))
}

#' @export
rev_cum_sum <- function(x) {
    l <- length(x)
    return(cumsum((x)[l:1])[l:1])
}

#' @export
compute_bounds <- function(x, sigma.ratio.lim) {
    l <- length(x)
    sel <- which(x >= sigma.ratio.lim[1] / (sigma.ratio.lim[1] + 1) * x[l])
    if (length(sel) > 0) {
        start <- max(1, min(sel))
    } else {
        start <- 1
    }
    sel <- which(x <= sigma.ratio.lim[2] / (sigma.ratio.lim[2] + 1) * x[l])
    if (length(sel) > 0) {
        end <- min(l - 1, max(sel))
    } else {
        end <- l - 1
    }
    return(list(start = start, end = end))
}

#' @export
bigauss.esti <- function(x, y, power = 1, do.plot = FALSE, truth = NA, sigma.ratio.lim = c(0.3, 3)) {
    sel <- which(y > 1e-10)
    if (length(sel) < 2) {
        to.return <- c(median(x), 1, 1, 0)
    } else {
        x <- x[sel]
        y <- y[sel]
        # 			sel<-order(x)
        # 			y<-y[sel]
        # 			x<-x[sel]

        y.0 <- y
        if (do.plot) plot(x, y)
        if (do.plot & !is.na(truth[1])) {
            true.y1 <- dnorm(x[x < truth[1]], mean = truth[1], sd = truth[2]) * truth[2] * truth[4]
            true.y2 <- dnorm(x[x >= truth[1]], mean = truth[1], sd = truth[3]) * truth[3] * truth[4]
            lines(x, c(true.y1, true.y2), col = "green")
        }
        max.y.0 <- max(y.0, na.rm = TRUE)
        y <- (y / max.y.0)^power

        l <- length(x)
        min.d <- min(diff(x))
        dx <- c(x[2] - x[1], (x[3:l] - x[1:(l - 2)]) / 2, x[l] - x[l - 1])
        if (l == 2) dx <- rep(diff(x), 2)
        dx[dx > 4 * min.d] <- 4 * min.d

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

        if (do.plot) abline(v = m)

        sel1 <- which(x < m)
        sel2 <- which(x >= m)
        s1 <- sqrt(sum((x[sel1] - m)^2 * y[sel1] * dx[sel1]) / sum(y[sel1] * dx[sel1]))
        s2 <- sqrt(sum((x[sel2] - m)^2 * y[sel2] * dx[sel2]) / sum(y[sel2] * dx[sel2]))


        if (power != 1) {
            s1 <- s1 * sqrt(power)
            s2 <- s2 * sqrt(power)
        }

        d1 <- dnorm(x[sel1], sd = s1, mean = m)
        d2 <- dnorm(x[sel2], sd = s2, mean = m)
        d <- c(d1 * s1, d2 * s2) # notice this "density" doesnt integrate to 1. Rather it integrates to (s1+s2)/2
        y <- y.0

        dy.ratio <- d^2 * log(y / d)
        dy.ratio[is.na(dy.ratio)] <- 0
        dy.ratio[dy.ratio == -Inf] <- 0
        dy.ratio[dy.ratio == Inf] <- 0


        scale <- exp(sum(dy.ratio) / sum(d^2))

        if (do.plot) {
            fit.1 <- d * scale
            lines(x[y > 0], fit.1, col = "red")
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

#' @export
bigauss.mix <- function(x, y, power = 1, do.plot = FALSE, sigma.ratio.lim = c(0.1, 10), bw = c(15, 30, 60), eliminate = .05, max.iter = 25, estim.method, BIC.factor = 2) {
    all.bw <- bw[order(bw)]

    x.0 <- x
    y.0 <- y

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
        this.smooth <- ksmooth(x.0, y.0, kernel = "normal", bandwidth = bw)
        turns <- find.turn.point(this.smooth$y)
        pks <- this.smooth$x[turns$pks]
        vlys <- c(-Inf, this.smooth$x[turns$vlys], Inf)

        smoother.pk.rec[[bw.n]] <- pks
        smoother.vly.rec[[bw.n]] <- vlys
        if (length(pks) != last.num.pks) {
            last.num.pks <- length(pks)
            l <- length(x)
            dx <- c(x[2] - x[1], (x[3:l] - x[1:(l - 2)]) / 2, x[l] - x[l - 1])
            if (l == 2) dx <- rep(diff(x), 2)

            # initiation
            m <- s1 <- s2 <- delta <- pks
            for (i in 1:length(m))
            {
                sel.1 <- which(x >= max(vlys[vlys < m[i]]) & x < m[i])
                s1[i] <- sqrt(sum((x[sel.1] - m[i])^2 * y[sel.1] * dx[sel.1]) / sum(y[sel.1] * dx[sel.1]))

                sel.2 <- which(x >= m[i] & x < min(vlys[vlys > m[i]]))
                s2[i] <- sqrt(sum((x[sel.2] - m[i])^2 * y[sel.2] * dx[sel.2]) / sum(y[sel.2] * dx[sel.2]))

                delta[i] <- (sum(y[sel.1] * dx[sel.1]) + sum(y[sel.2] * dx[sel.2])) / ((sum(dnorm(x[sel.1], mean = m[i], sd = s1[i])) * s1[i] / 2) + (sum(dnorm(x[sel.2], mean = m[i], sd = s2[i])) * s2[i] / 2))
            }
            delta[is.na(delta)] <- 1e-10
            s1[is.na(s1)] <- 1e-10
            s2[is.na(s2)] <- 1e-10


            fit <- matrix(0, ncol = length(m), nrow = length(x)) # this is the matrix of fitted values

            this.change <- Inf
            counter <- 0

            while (this.change > 0.1 & counter <= max.iter) {
                counter <- counter + 1
                old.m <- m

                # E step
                cuts <- c(-Inf, m, Inf)
                for (j in 2:length(cuts))
                {
                    sel <- which(x >= cuts[j - 1] & x < cuts[j])
                    use.s1 <- which(1:length(m) >= (j - 1))
                    s.to.use <- s2
                    s.to.use[use.s1] <- s1[use.s1]
                    for (i in 1:ncol(fit))
                    {
                        fit[sel, i] <- dnorm(x[sel], mean = m[i], sd = s.to.use[i]) * s.to.use[i] * delta[i]
                    }
                }
                fit[is.na(fit)] <- 0
                sum.fit <- apply(fit, 1, sum)

                # Elimination step
                fit <- fit / sum.fit
                fit2 <- fit * y
                perc.explained <- apply(fit2, 2, sum) / sum(y)
                max.erase <- max(1, round(length(perc.explained) / 5))

                to.erase <- which(perc.explained <= min(eliminate, perc.explained[order(perc.explained, na.last = FALSE)[max.erase]]))


                if (length(to.erase) > 0) {
                    m <- m[-to.erase]
                    s1 <- s1[-to.erase]
                    s2 <- s2[-to.erase]
                    delta <- delta[-to.erase]
                    fit <- fit[, -to.erase]
                    if (is.null(ncol(fit))) fit <- matrix(fit, ncol = 1)
                    sum.fit <- apply(fit, 1, sum)
                    fit <- fit / sum.fit
                    old.m <- old.m[-to.erase]
                }

                # M step
                for (i in 1:length(m))
                {
                    this.y <- y * fit[, i]
                    if (estim.method == "moment") {
                        this.fit <- bigauss.esti(x, this.y, power = power, do.plot = FALSE, sigma.ratio.lim = sigma.ratio.lim)
                    } else {
                        this.fit <- bigauss.esti.EM(x, this.y, power = power, do.plot = FALSE, sigma.ratio.lim = sigma.ratio.lim)
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
                sel <- which(x >= cuts[j - 1] & x < cuts[j])
                use.s1 <- which(1:length(m) >= (j - 1))
                s.to.use <- s2
                s.to.use[use.s1] <- s1[use.s1]
                for (i in 1:ncol(fit))
                {
                    if (s.to.use[i] != 0) {
                        fit[sel, i] <- dnorm(x[sel], mean = m[i], sd = s.to.use[i]) * s.to.use[i] * delta[i]
                    }
                }
            }

            if (do.plot) {
                plot(x, y, cex = .1, main = paste("bw=", bw))
                sum.fit <- apply(fit, 1, sum)
                lines(x, sum.fit)
                abline(v = m)
                cols <- c("red", "green", "blue", "cyan", "brown", "black", rep("grey", 100))
                for (i in 1:length(m))
                {
                    lines(x, fit[, i], col = cols[i])
                }
            }
            area <- delta * (s1 + s2) / 2
            rss <- sum((y - apply(fit, 1, sum))^2)
            l <- length(x)
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
    if (length(sel) > 1) sel <- sel[which(all.bw[sel] == max(all.bw[sel]))]
    rec <- new("list")
    rec$param <- results[[sel]]
    rec$smoother.pks <- smoother.pk.rec
    rec$smoother.vlys <- smoother.vly.rec
    rec$all.param <- results
    rec$bic <- bic.rec
    return(rec)
}

#' @export
normix <- function(that.curve, pks, vlys, ignore = 0.1, max.iter = 50, prob.cut = 1e-2) {
    x <- that.curve[, 1]
    y <- that.curve[, 2]

    if (length(pks) == 1) {
        miu <- sum(x * y) / sum(y)
        sigma <- sqrt(sum(y * (x - miu)^2) / sum(y))
        fitted <- dnorm(x, mean = miu, sd = sigma)
        this.sel <- y > 0 & fitted / dnorm(miu, mean = miu, sd = sigma) > prob.cut
        sc <- exp(sum(fitted[this.sel]^2 * log(y[this.sel] / fitted[this.sel]) / sum(fitted[this.sel]^2)))
    } else {
        pks <- pks[order(pks)]
        vlys <- vlys[order(vlys)]
        l <- length(pks)
        miu <- sigma <- sc <- pks
        w <- matrix(0, nrow = l, ncol = length(x))

        for (m in 1:l)
        {
            this.low <- max(vlys[vlys <= pks[m]])
            this.high <- min(vlys[vlys >= pks[m]])
            this.x <- x[x >= this.low & x <= this.high]
            this.y <- y[x >= this.low & x <= this.high]

            miu[m] <- sum(this.x * this.y) / sum(this.y)
            # if(sum(this.y>0) > 1)
            # {
            sigma[m] <- sqrt(sum(this.y * (this.x - miu[m])^2) / sum(this.y))
            # }else{
            #  sigma[m]<-mean(diff(this.x))/2
            # }
            fitted <- dnorm(this.x, mean = miu[m], sd = sigma[m])
            this.sel <- this.y > 0 & fitted / dnorm(miu[m], mean = miu[m], sd = sigma[m]) > prob.cut
            sc[m] <- exp(sum(fitted[this.sel]^2 * log(this.y[this.sel] / fitted[this.sel]) / sum(fitted[this.sel]^2)))
            # sc[m]<-lm(this.y[this.sel]~fitted[this.sel]+0)$coef
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
                miu <- sum(x * y) / sum(y)
                sigma <- sqrt(sum(y * (x - miu)^2) / sum(y))
                fitted <- dnorm(x, mean = miu, sd = sigma)
                this.sel <- y > 0 & fitted / dnorm(miu, mean = miu, sd = sigma) > prob.cut
                sc <- exp(sum(fitted[this.sel]^2 * log(y[this.sel] / fitted[this.sel]) / sum(fitted[this.sel]^2)))
                # sc<-lm(y[this.sel]~fitted[this.sel]+0)$coef
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

            if (sum(is.na(w)) > 0) break

            for (m in 1:l)
            {
                this.y <- y * w[m, ]
                miu[m] <- sum(x * this.y) / sum(this.y)
                sigma[m] <- sqrt(sum(this.y * (x - miu[m])^2) / sum(this.y))
                fitted <- dnorm(x, mean = miu[m], sd = sigma[m])
                this.sel <- this.y > 0 & fitted / dnorm(miu[m], mean = miu[m], sd = sigma[m]) > prob.cut
                sc[m] <- exp(sum(fitted[this.sel]^2 * log(this.y[this.sel] / fitted[this.sel]) / sum(fitted[this.sel]^2)))
                # sc[m]<-lm(this.y[this.sel]~fitted[this.sel]+0)$coef
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

#' @export
normix.bic <- function(x, y, power = 2, do.plot = FALSE, bw = c(15, 30, 60), eliminate = .05, max.iter = 50, BIC.factor = 2) {
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
            aaa <- normix(cbind(x, y), pks = pks, vlys = vlys, ignore = eliminate, max.iter = max.iter)

            total.fit <- x * 0
            for (i in 1:nrow(aaa))
            {
                total.fit <- total.fit + dnorm(x, mean = aaa[i, 1], sd = aaa[i, 2]) * aaa[i, 3]
            }

            if (do.plot) {
                plot(x, y, cex = .1, main = paste("bw=", bw))
                abline(v = aaa[, 1])
                cols <- c("red", "green", "blue", "cyan", "brown", "black", rep("grey", 100))
                for (i in 1:nrow(aaa))
                {
                    lines(x, dnorm(x, mean = aaa[i, 1], sd = aaa[i, 2]) * aaa[i, 3], col = cols[i])
                }
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
    if (length(sel) > 1) sel <- sel[which(all.bw[sel] == max(all.bw[sel]))]
    rec <- new("list")
    rec$param <- results[[sel]]
    rec$smoother.pks <- smoother.pk.rec
    rec$smoother.vlys <- smoother.vly.rec
    rec$all.param <- results
    rec$bic <- bic.rec
    return(rec)
}

#' Generate feature table from noise-removed LC/MS profile
#' 
#' @description
#' Each LC/MS profile is first processed by the function proc.cdf() to remove noise and reduce data size. A matrix containing m/z 
#' value, retention time, intensity, and group number is output from proc.cdf(). This matrix is then fed to the function 
#' prof.to.features() to generate a feature list. Every detected feature is summarized into a single row in the output matrix from this function.
#' 
#' @param feature_table The matrix output from proc.cdf(). It contains columns of m/z value, retention time, intensity and group number.
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
prof.to.features <- function(feature_table,
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

    feature_table <- preprocess_feature_table(feature_table)

    bws <- preprocess_bandwidth(min.bw, max.bw, feature_table)
    min.bw <- bws[["min.bw"]]
    max.bw <- bws[["max.bw"]]

    base.curve <- compute_base_curve(feature_table[, "rt"])
    all.times <- compute_all_times(base.curve)

    keys <- c("mz", "pos", "sd1", "sd2", "area")
    this.features <- matrix(0, nrow = 1, ncol = length(keys), dimnames = list(NULL, keys))

    feature_groups <- split(feature_table, feature_table$group_number)
    for (i in seq_along(feature_groups))
    {
        this <- feature_groups[[i]]
        this <- this[order(this[, "rt"]), ]

        mz.sd.rec <- c(NA, sd(this[, "mz"]))

        nrow_this <- nrow(this)
        if (between(nrow_this, 2, 10)) {
            this.inte <- interpol.area(this[, "rt"], this[, "intensity"], base.curve[, 1], all.times)
            chr_peak_shape <- c(median(this[, "mz"]), median(this[, "rt"]), sd(this[, "rt"]), sd(this[, "rt"]), this.inte)
            this.features <- rbind(this.features, chr_peak_shape)
        }

        if (nrow_this < 2) {
            this.time.weights <- all.times[which(base.curve[, 1] %in% this[2])]
            chr_peak_shape <- c(this[1], this[2], NA, NA, this[3] * this.time.weights)
            this.features <- rbind(this.features, chr_peak_shape)
        }

        if (nrow_this > 10) {
            this.span <- range(this[, "rt"])
            this.curve <- base.curve[base.curve[, "base_curve"] >= this.span[1] & base.curve[, "base_curve"] <= this.span[2], ]
            this.curve[this.curve[, "base_curve"] %in% this[, "rt"], 2] <- this[, "intensity"]

            bw <- min(max(bandwidth * (this.span[2] - this.span[1]), min.bw), max.bw)
            bw <- seq(bw, 2 * bw, length.out = 3)

            if (bw[1] > 1.5 * min.bw) {
                bw <- c(max(min.bw, bw[1] / 2), bw)
            }

            if (shape.model == "Gaussian") {
                chr_peak_shape <- compute_gaussian_peak_shape(this.curve, power, bw, component.eliminate, BIC.factor)
            } else {
                chr_peak_shape <- bigauss.mix(this.curve[, "base_curve"], this.curve[, 2], sigma.ratio.lim = sigma.ratio.lim, bw = bw, power = power, estim.method = estim.method, eliminate = component.eliminate, BIC.factor = BIC.factor)$param[, c(1, 2, 3, 5)]
            }

            if (is.null(nrow(chr_peak_shape))) {
                this.features <- rbind(this.features, c(median(this[, "mz"]), chr_peak_shape))
            } else {
                for (m in 1:nrow(chr_peak_shape))
                {
                    this.d <- abs(this[, "rt"] - chr_peak_shape[m, 1])
                    this.features <- rbind(this.features, c(mean(this[which(this.d == min(this.d)), 1]), chr_peak_shape[m, ]))
                }
            }
        }
    }
    this.features <- this.features[-1, ]
    this.features <- this.features[order(this.features[, "mz"], this.features[, "pos"]), ]
    this.features <- this.features[which(apply(this.features[, c("sd1", "sd2")], 1, min) > sd.cut[1] & apply(this.features[, c("sd1", "sd2")], 1, max) < sd.cut[2]), ]
    rownames(this.features) <- NULL

    if (do.plot) {
        par(mfrow = c(2, 2))
        plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
        text(x = 0, y = 0, "Estimate peak \n area/location", cex = 1.5)
        hist(mz.sd.rec, xlab = "m/z SD", ylab = "Frequency", main = "m/z SD distribution")
        hist(c(this.features[, "sd1"], this.features[, "sd2"]), xlab = "Retention time SD", ylab = "Frequency", main = "Retention time SD distribution")
        hist(log10(this.features[, "area"]), xlab = "peak strength (log scale)", ylab = "Frequency", main = "Peak strength distribution")
    }

    return(this.features)
}
