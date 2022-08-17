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
#' @examples
#' rm.ridge(rt, intensi, bw)
rm.ridge <-
function(x,y2, bw)
{
    sel<-which(y2<quantile(y2, 0.75))
    max.x.sel<-max(x[sel])
    min.x.sel<-min(x[sel])
    
    in.sel<-which(between(x, min.x.sel, max.x.sel))
    over.sel<-which(x>max.x.sel)
    under.sel<-which(x<min.x.sel)
    
    
    this.s<-ksmooth(x[sel],y2[sel],x.points=x[in.sel],kernel="normal",bandwidth=bw)
    if(sum(is.na(this.s$y))>0) return(y2)
    
    y2[in.sel]<-y2[in.sel]-this.s$y
    y2[over.sel]<-y2[over.sel]-this.s$y[which(this.s$x==max(this.s$x))[1]]
    y2[under.sel]<-y2[under.sel]-this.s$y[which(this.s$x==min(this.s$x))[1]]
    
    y2[y2<0]<-0
    return(y2)
}
