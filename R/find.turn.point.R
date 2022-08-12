#' @description
#' Compute local maxima turn points. 
#' @param y The y values of a curve in x-y plane.
#' @param ties.method specifies the method rank uses to break ties.
#' @return boolean row with local maxima turn point.
#' @export
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

    return(v)
}

#' @description
#' Compute maxima and minima turn points. 
#' @param y The y values of a curve in x-y plane.
#' @return boolean row with local maxima and minima turn points.
#' @export
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
find.turn.point <- function(y) {
    y <- y[!is.na(y)] # filter NA values

    if (length(unique(y)) == 1) { # if exactly one distinct value
        middle_index <- round(length(y) / 2) # get mid index
        start_and_end <- c(1, length(y)) # get first and last index
        return(list(pks = middle_index, vlys = start_and_end))
    } else {
        b <- msExtrema(y)

        pks <- which(b$index.max)
        vlys <- which(b$index.min)

        if (pks[1] != 1) {
            vlys <- c(1, vlys)
        }

        if (pks[length(pks)] != length(y)) {
            vlys <- c(vlys, length(y))
        }

        if (length(pks) == 1) {
            vlys <- c(1, length(y))
        }

        x <- new("list")
        x$pks <- pks
        x$vlys <- vlys
        return(x)
    }
}
