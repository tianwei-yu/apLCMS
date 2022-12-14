#' Removing long ridges at the same m/z.
#' 
#' @description
#' This is an internal function. It substracts a background estimated through kernel smoothing when an EIC continuously 
#' span more than half the retention time range.
#' @param x Retention time vector.
#' @param y2 Intensity vector.
#' @param bw Bandwidth for the kernel smoother. A very wide one is used here.
#' @return A vector of intensity value is returned.
#' @importFrom dplyr between
#' @export
rm.ridge <- function(x,y2, bw) {
    
    sel <- which(y2<quantile(y2, 0.75))
    sel_max_min <- tibble::tibble(max_sel = max(x[sel]), min_sel = min(x[sel]))
    max.x.sel <- sel_max_min$max_sel
    min.x.sel <- sel_max_min$min_sel
     
    in.sel <- which(between(x, min.x.sel, max.x.sel))
    sel_over_under <- tibble::tibble(over_sel = which(x > max.x.sel), under_sel = which(x < min.x.sel))
    over.sel <- sel_over_under$over_sel 
    under.sel <- sel_over_under$under_sel 
    
    
    this.s <- ksmooth(x[sel], y2[sel], x.points = x[in.sel], kernel = "normal", bandwidth = bw)
    if(sum(is.na(this.s$y)) > 0) return(y2)
    
    y2[in.sel] <- y2[in.sel] - this.s$y
    y2[over.sel] <- y2[over.sel]-this.s$y[which(this.s$x==max(this.s$x))[1]]
    y2[under.sel] <- y2[under.sel]-this.s$y[which(this.s$x==min(this.s$x))[1]]
    
    y2[y2<0] <- 0
    return(y2)
}
