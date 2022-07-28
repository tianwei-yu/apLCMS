#' Compute matches between mz array and specific mass value with a tolerance.
#' @param x The mz array for which to compute the matching.
#' @param known.mz The mz value with which to match.
#' @param match.tol.ppm Matching tolerance in ppm.
#' @return  Binary vector, 1 indicating a match, 0 a mismatch.
#' @examples
#' mass.match(aligned.ftrs, known.table, match.tol.ppm)
mass.match <-
function(x, known.mz, match.tol.ppm=5)
{
    mass.matched.pos<-rep(0, length(x))
    for(i in 1:length(x))
    {
        this.d<-abs((x[i]-known.mz)/x[i])
        if(min(this.d) < match.tol.ppm/1e6) mass.matched.pos[i]<-1
    }
    return(mass.matched.pos)
}
