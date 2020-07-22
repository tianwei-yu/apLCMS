\name{prof.learn}
\alias{prof.learn}
\docType{data}
\title{ Sample profile data after noise filtration by the machine learning approach. }
\description{
  A list object containing 4 matrices. Each matrix is from an LC/MS profile. 
}
\usage{data(prof)}
\format{
Each matrix contains 4 columns: m/z, retention time, intensity, and group number. 
}
\details{}
\source{Data from Dean Jones lab, Emory University School of Medicine.}
\references{}
\examples{
data(prof.learn)
this.feature<-prof.to.features(prof.learn[[1]][[1]])
}
\keyword{datasets}
