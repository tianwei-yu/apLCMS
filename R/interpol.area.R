#' Interpolate missing intensities and calculate the area for a single EIC.
#' 
#' This is an internal function.
#' 
#' @param x the positions of x(retention time) where non-NA y is observed.
#' @param y the observed intensities.
#' @param all.x all possible x(retention time) in the LCMS profile.
#' @param all.w the "footprint" of each measured retention time, used as weight for the corresponding y.
#' @return The area is returned.
#' @export
#' @examples
#' interpol.area(x, y, all.x, all.w)
interpol.area <-
function(x, y, all.x, all.w) # x is retention time, all.w is weight of each all.x, y is observed intensities
{
    r<-range(x)
    indices <- which(between(all.x, r[1], r[2]))
    
    x2<-all.x[indices]
    w<-all.w[indices]
    
    y2<-approx(x, y, xout=x2, method="linear")$y
    return(sum(w*y2))
}
