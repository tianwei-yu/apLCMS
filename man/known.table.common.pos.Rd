\name{known.table.common.pos}
\alias{known.table.common.pos}
\docType{data}
\title{
A known feature table based on HMDB.
}
\description{
This table was compiled from HMDB metabolites. It contains only the H+, Na+, K+ and NH4+ derivative of known metabolites. It suites the need for analyzing LC/MS data generated with anion exchange chromatography. 
}
\usage{data(known.table.common.pos)}
\format{
A data frame containing the known metabolite ions and previously found features. It contains 18 columns: 
	"chemical_formula": the chemical formula if knonw; 
	"HMDB_ID": HMDB ID if known; 
	"KEGG_compound_ID": KEGG compound ID if known; 
	"neutral.mass": the neutral mass if known: 
	"ion.type": the ion form, such as H+, Na+, ..., if known; 
	"m.z": m/z value, either theoretical for known metabolites, or mean observed value for unknown but previously found features; 
	"Number_profiles_processed": the total number of LC/MS profiles that were used to build this database;
	"Percent_found": in what percentage was this feature found historically amount all data processed in building this database;
	"mz_min": the minimum m/z value observed for this feature;
	"mz_max": the maximum m/z value observed for this feature;
	"RT_mean": the mean retention time observed for this feature;
	"RT_sd": the standard deviation of retention time observed for this feature;
	"RT_min": the minimum retention time observed for this feature;
	"RT_max": the maximum retention time observed for this feature;
	"int_mean.log.": the mean log intensity observed for this feature; 
	"int_sd.log.": the standard deviation of log intensity observed for this feature;
	"int_min.log.": the minimum log intensity observed for this feature;
	"int_max.log.": the maximum log intensity observed for this feature;
}
\details{}
\source{
Wishart, D. S., et al. (2009). HMDB: a knowledgebase for the human metabolome. Nucleic Acids Res 37, D603-10.}
\references{}
\examples{
data(known.table.common.pos)
}
\keyword{datasets}