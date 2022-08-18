#' Internal function: Updates the information of a feature for the known feature table.
#' 
#' @description
#' The function takes the information about the feature in the known feature table (if available), and updates it using the 
#' information found in the current dataset.
#' @param existing.row The existing row in the known feature table.
#' @param ftrs.row The row of the matched feature in the new aligned feature table.
#' @param rt.row The row of the matched feature in the new retention time table of aligned features.
#' @return A vector, the updated row for the known feature table.
#' @examples
#' peak.characterize(existing.row = NA, ftrs.row, rt.row)
peak.characterize <-
function(existing.row=NA, ftrs.row, rt.row)
{
    merge.new<-function(mean0, sd0, min0, max0, n, x)
    {
        x<-x[!is.na(x)]
        if(n<=1)
        {
            if(n==1) x<-c(x,mean0)
            mean1<-mean(x)
            sd1<-sd(x)
            min1<-min(x, Inf)
            max1<-max(x, -Inf)
        }else{
            m<-length(x)
            mean1<-sum(mean0*n, x)/(n+m)
            sd1<-sqrt((n*(mean0-mean1)^2+sum((x-mean1)^2)+(n-1)*sd0^2)/(m+n-1))
            min1<-min(min0, x)
            max1<-max(max0, x)
        }
        return(c(mean1, sd1, min1, max1))
    }
    
    
    ftrs.row[5:length(ftrs.row)]<-log10(ftrs.row[5:length(ftrs.row)]+1)
    ftrs.row[ftrs.row==0]<-NA
    if(length(existing.row)==1)
    {
        existing.row<-rep(NA, 18)
        existing.row[6]<-ftrs.row[1]
    }
    
    n<-round(as.numeric(existing.row[7])*as.numeric(existing.row[8])) #times found in previous experiments
    if(is.na(n)) n<-0
    m<-sum(!is.na(rt.row[5:length(rt.row)])) #times found in current experiment
    
    existing.row[7]<-sum(as.numeric(existing.row[7]), length(ftrs.row)-4, na.rm=T)
    existing.row[8]<-(n+m)/as.numeric(existing.row[7])
    existing.row[9]<-min(as.numeric(existing.row[6]), as.numeric(existing.row[9]), ftrs.row[3],na.rm=T)
    existing.row[10]<-max(as.numeric(existing.row[6]), as.numeric(existing.row[10]), ftrs.row[4],na.rm=T)
    
    this<-merge.new(as.numeric(existing.row[11]),as.numeric(existing.row[12]),as.numeric(existing.row[13]),as.numeric(existing.row[14]), n, rt.row[5:length(rt.row)])
    existing.row[11:14]<-this
    
    this<-merge.new(as.numeric(existing.row[15]),as.numeric(existing.row[16]),as.numeric(existing.row[17]),as.numeric(existing.row[18]), n, ftrs.row[5:length(ftrs.row)])
    existing.row[15:18]<-this
    
    return(existing.row)
}
