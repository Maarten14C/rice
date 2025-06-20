% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sets.R
\name{pool}
\alias{pool}
\title{Test if a set of radiocarbon dates can be combined}
\usage{
pool(
  y,
  er,
  deltaR = 0,
  deltaSTD = 0,
  threshold = 0.05,
  roundby = 1,
  talk = TRUE
)
}
\arguments{
\item{y}{The set of radiocarbon dates to be tested}

\item{er}{The lab errors of the radiocarbon dates}

\item{deltaR}{Age offset (e.g. for marine samples).}

\item{deltaSTD}{Uncertainty of the age offset (1 standard deviation).}

\item{threshold}{Probability threshold above which chisquare values are considered acceptable (between 0 and 1; default \code{threshold=0.05}).}

\item{roundby}{Rounding of the reported mean, chisquare and and p-value. Defaults to \code{roundby=1}.}

\item{talk}{It's better than staying silent.}
}
\value{
The pooled mean and error if the p-value is above the threshold - a warning if it is not.
}
\description{
Calculate the (chi-square) probability that a set of radiocarbon dates is consistent, i.e. that it can be assumed that they all pertain to the same true radiocarbon age (and thus to the same calendar age - note though that sometimes multiple calendar ages obtain the same C14 age). The function calculates the differences (chi2 value) and finds the corresponding p-value. If the chi2 values is sufficiently small, then the p-value is sufficiently large (above the threshold), and the pooled mean is calculated and returned. If the scatter is too large, no pooled mean is calculated.
}
\details{
This follows the calculations of Ward and Wilson (1978; Archaeometry 20: 19-31 <doi:10.1111/j.1475-4754.1978.tb00208.x>) and should only be used for multiple dates that stem from the same sample (e.g., multiple measurements on a single bone). It cannot be used to test if multiple dates from multiple samples pertain to the same event. Since the assumption is that all measurements stem from the same event, we can assume that they all share the same C14 age (since any calBP age will have an associated IntCal C14 age).
}
\examples{
  data(shroud)
  pool(shroud$y,shroud$er)
  Zu <- grep("ETH", shroud$ID) # Zurich lab only
  pool(shroud$y[Zu],shroud$er[Zu])
}
\author{
Maarten Blaauw
}
