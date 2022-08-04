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
combine.seq.3 <- function(features) {
    l <- nrow(features)
    breaks <- c(0, which(features$labels[1:(l - 1)] != features$labels[2:l]), l)
    new_table <- data.frame(masses = rep(0, length(breaks) - 1), labels = unique(features$labels), intensi = rep(0, length(breaks) - 1))
    
    for (i in 1:(length(breaks) - 1)) {
        start <- breaks[i] + 1
        end <- breaks[i + 1]
        mz <- features$masses[start:end]
        ints <- features$intensi[start:end]

        new_table$intensi[i] <- sum(ints)
        new_table$masses[i] <- median(mz[which.max(ints)])
    }

    return(new_table)
}

combine.seq.3_old <-
    function(a, mz, inte) ### the input need to be pre-ordered by a
    {
        l <- length(a)
        breaks <- c(0, which(a[1:(l - 1)] != a[2:l]), l)
        new.int <- new.mz <- rep(0, length(breaks) - 1)

        for (i in 1:(length(breaks) - 1)) {
            this.int <- inte[(breaks[i] + 1):breaks[i + 1]]
            this.mz <- mz[(breaks[i] + 1):breaks[i + 1]]
            new.int[i] <- sum(this.int)
            new.mz[i] <- median(this.mz[which(this.int == max(this.int))])
        }
        new.a <- unique(a)
        return(cbind(new.mz, new.a, new.int))
    }
