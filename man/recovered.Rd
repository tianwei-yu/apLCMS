\name{recovered}
\alias{recovered}
\docType{data}
\title{Sample data after weak signal recovery}
\description{This is a list object with four items, each from a sample LC/MS profile.}
\usage{data(recovered)}
\format{
  The format is:
List of 4, each item is in turn a list of 4.
	$this.ftrs		: the signal strength of each feature identified.
	$this.times		: the retention time of each feature identified.
	$this.f1		: A matrix with m/z, retention time, retention time spread, and peak area
	this.f2			: A matrix with m/z, retention time, retention time spread, peak area, file identifier, and feature identifier
}
\details{}
\source{Dean Jones lab. Emory University, School of Medicine.}
\references{}
\examples{
data(recovered)
}
\keyword{datasets}