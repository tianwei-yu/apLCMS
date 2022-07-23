find.turn.point <- function(values) {
  
  find_extrema <- function(values, ties.method) {
    # this might not work correctly
    z <- embed(rev(as.vector(c(-Inf, values,-Inf))), dim = 3)
    z <- z[rev(seq(nrow(z))),]
    v <- max.col(z, ties.method = ties.method) == 2
    return(v)
  }
  
  analyse_extrema <- function(values) {
    minima <- find_extrema(values, ties.method = "first") # find peaks
    maxima <- find_extrema(-values, ties.method = "last") # find valleys
    # remove overlaps
    maxima.indices <- minima & !maxima
    minima.indices <- maxima & !minima
    return(list(minima.indices = minima.indices, maxima.indices = maxima.indices))
  }
  
  values <- values[!is.na(values)]
  # flat curve => peak is middle and valleys are the first and the last
  if (length(unique(values)) == 1) {
    x <- new("list")
    x$pks <- round(length(values) / 2)
    x$vlys <- c(1, length(values))
    return(x)
  }
  
  extrema <- analyse_extrema(values)
  maxima <- which(extrema$maxima.indices)
  minima <- which(extrema$minima.indices)
  # if the first or the last index are not maxima, added them to minima
  if (maxima[1] != 1)
    minima <- c(1, minima)
  if (maxima[length(maxima)] != length(values))
    minima <- c(minima, length(values))
  # if there is only one maximum, the minima are the first and the last index
  # (probably because Gaussian was used)
  if (length(maxima) == 1)
    minima <- c(1, length(values))
  x <- new("list")
  x$pks <- maxima
  x$vlys <- minima
  return(x)
}
