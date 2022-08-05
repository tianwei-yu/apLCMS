#' Find peaks and valleys of a curve.
#'
#' @description
#' This is an internal function which finds the peaks and valleys of a smooth curve.
#'
#' @param y The y values of a curve in x-y plane.
#' @return A list object:
#' \itemize{
#'   \item pks - The peak positions
#'   \item vlys - The valley positions
#' }
#' @export
#' @examples
#' find.turn.point(y)
library(pastecs)
find.turn.point <- function(y) {
    y <- y[!is.na(y)]                             # filter NA values
    if (length(unique(y)) == 1) {                 # if exactly one distinct value
        middle_index <- round(length(y) / 2)      # get mid index
        start_and_end <- c(1, length(y))          # get first and last index
        return(list(pks = middle_index, vlys = start_and_end))
    } else {
        list_tp <- turnpoints(y)        
        pks <- which(list_tp$peaks)  
        vlys <- which(list_tp$pits)
        if (length(pks) == 1) {   
            vlys <- c(1, list_tp$n)
        }
        x <- new("list")
        x$pks <- pks
        x$vlys <- vlys
    return(x)
    }
}
