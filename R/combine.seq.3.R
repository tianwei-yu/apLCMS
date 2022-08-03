#' An internal function.
#' 
#' This is a internal function.
#' 
#' @param table dataframe of retention time, m/z ratio, signal strength.
#' @return returns
#' \itemize{
#'   \item masses - m/z ratio
#'   \item labels - retention time
#'   \item intensi - signal strength
#' }
#' @export
#' @examples
#' combine.seq.3(table)
combine.seq.3 <- function(table) {
    l <- nrow(table)
    breaks <- c(0, which(table$labels[1:(l - 1)] != table$labels[2:l]), l)
    new_table <- data.frame(masses = rep(0, length(breaks) - 1), labels = unique(table$labels), intensi = rep(0, length(breaks) - 1))
    
    for (i in 1:(length(breaks) - 1)) {
        start <- breaks[i] + 1
        end <- breaks[i + 1]

        this_table <- data.frame(intensi = table$intensi[start:end], masses = table$masses[start:end])

        new_table$intensi[i] <- sum(this_table$intensi)
        new_table$masses[i] <- median(this_table$masses[which(this_table$intensi == max(this_table$intensi))])
    }

    return(new_table)
}