\name{unsupervised}
\alias{unsupervised}
\title{Runs features extraction in unsupervised mode.}
\description{TBA}
\usage{
unsupervised(filenames,
  min_exp = 2,
  min_pres = 0.5,
  min_run = 12,
  mz_tol = 1e-05,
  baseline_correct = 0,
  baseline_correct_noise_percentile = 0.05,
  shape_model = "bi-Gaussian",
  BIC_factor = 2,
  peak_estim_method = "moment",
  min_bandwidth = NA,
  max_bandwidth = NA,
  sd_cut = c(0.01, 500),
  sigma_ratio_lim = c(0.01, 100),
  component_eliminate = 0.01,
  moment_power = 1,
  align_mz_tol = NA,
  align_chr_tol = NA,
  max_align_mz_diff = 0.01,
  recover_mz_range = NA,
  recover_chr_range = NA,
  use_observed_range = TRUE,
  recover_min_count = 3,
  intensity_weighted = FALSE,
  cluster = 4
)
}
\arguments{
  \item{filenames}{The CDF file names.}
  \item{min_exp}{TBA}
  \item{min_pres}{TBA}
  \item{min_run}{Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.}
  \item{mz_tol}{m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.}
  \item{baseline_correct}{After grouping the observations, the highest intensity in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the height of the noise groups.}
  \item{baseline_correct_noise_percentile}{The perenctile of signal strength of those EIC that don't pass the run filter, to be used as the baseline threshold of signal strength.}
  \item{shape_model}{The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian". When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".}
  \item{BIC_factor}{The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1, models with more peaks are penalized more.}
  \item{peak_estim_method}{The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.}
  \item{min_bandwidth}{The minimum bandwidth to use in the kernel smoother.}
  \item{max_bandwidth}{The maximum bandwidth to use in the kernel smoother.}
  \item{sd_cut}{A vector of two. Features with standard deviation outside the range defined by the two numbers are eliminated.}
  \item{sigma_ratio_lim}{A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation and the righ-standard deviation of the bi-Gaussian fuction used to fit the data.}
  \item{component_eliminate}{In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a proportion of intensities less than this value, the component will be ignored.}
  \item{moment_power}{The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.}
  \item{align_mz_tol}{The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.}
  \item{align_chr_tol}{The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for the tolerance level based on the data.}
  \item{max_align_mz_diff}{As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.}
  \item{recover_mz_range}{The m/z around the feature m/z to search for observations. The default value is NA, in which case 1.5 times the m/z tolerance in the aligned object will be used.}
  \item{recover_chr_range}{The retention time around the feature retention time to search for observations. The default value is NA, in which case 0.5 times the retention time tolerance in the aligned object will be used.}
  \item{use_observed_range}{If the value is TRUE, the actual range of the observed locations of the feature in all the spectra will be used.}
  \item{recover_min_count}{Minimum number of raw data points to support a recovery.}
  \item{intensity_weighted}{TBA}
  \item{cluster}{Number of requered processor cores.}
}
\value{ TBA }
\author{ Jiri Novotny <xtracko@fi.muni.cz> }
\keyword{ models }
