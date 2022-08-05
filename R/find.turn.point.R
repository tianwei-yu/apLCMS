find_local_maxima <- function(y, ties.method) {
    padded_y <- rev(as.vector(c(-Inf, y, -Inf)))

    # each row is 3 consecutive values in descending order
    # rows are sorted in ascending order
    z <- embed(padded_y, dim = 3)
    
    # reverse the row ordering
    # first column is equal y
    z <- z[rev(seq(nrow(z))), ]

    # row where the max is in the middle is a turn point
    v <- max.col(z, ties.method = ties.method) == 2

    # alt <- diff(sign(diff(y))) == -2
    return(v)
}

msExtrema <- function(y) {
    index1 <- find_local_maxima(y, ties.method = "first")
    index2 <- find_local_maxima(-y, ties.method = "last")

    # this is some sort of safety mechanism to protect against numerical issues
    index.max <- index1 & !index2
    index.min <- index2 & !index1

    list(index.max = index.max, index.min = index.min)
}

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
    
    y <- y[!is.na(y)] # filter NA values
    m <- turnpoints(y)  # It gives the index.max peak

    if (length(unique(y)) == 1) { # if exactly one distinct value
        middle_index <- round(length(y) / 2) # get mid index
        start_and_end <- c(1, length(y)) # get first and last index
        return(list(pks = middle_index, vlys = start_and_end))
    } else {
        b <- msExtrema(y)  # Boolean list, TRUE for the extreme and peak values

        pks <- which(m$peaks) #which(b$index.max)    # gives 1 index -> 258
        #vlys <- which(b$index.min)   # gives 2 indices -> 1, 512
        
        # These If redefine the vlys, not sure it is correct !
        # if (pks[1] != 1) {    # this is TRUE 
        #     vlys <- c(1, vlys)
        # }

        # if (pks[length(pks)] != length(y)) {  # This is TRUE
        #     vlys <- c(vlys, length(y))
        # }

        if (length(pks) == 1) {   #This is TRUE -> final vlys values are assigned here 
            vlys <- c(1, m$n)#c(1, length(y))
        }

        x <- new("list")
        x$pks <- pks
        x$vlys <- vlys
    return(x)
    }
}
