% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realms.R
\name{b2ktoBCAD}
\alias{b2ktoBCAD}
\title{calculate BC/AD ages from b2k ages}
\usage{
b2ktoBCAD(x, zero = TRUE)
}
\arguments{
\item{x}{The b2k age(s) to be translated into BC/AD ages.}

\item{zero}{Whether or not zero BC/AD should be included. Defaults to \code{zero=TRUE}.}
}
\value{
The BC/AD age(s). BC ages are negative, AD ages are positive.
}
\description{
calculate BC/AD ages from b2k ages
}
\details{
Turn b2k ages (popular in the ice core community) into BC/AD (or BCE/CE). Negative ages indicate BC, positive ages AD. Since the Gregorian and Julian calendars do not include 0 BCAD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from -1 (i.e., 1 BC) to 1 AD. Other calendars, such as the astronomical one, do include zero. The often-used BCE/CE ages are equivalent to BC/AD.
}
\examples{
 b2ktoBCAD(0)
 b2ktoBCAD(1990:2010, zero=TRUE)
 b2ktoBCAD(1990:2010, zero=FALSE)
}
