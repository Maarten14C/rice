% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marine.R
\name{find.shells}
\alias{find.shells}
\title{Find nearby shell-derived dR values}
\usage{
find.shells(
  longitude,
  latitude,
  nearest = 50,
  colour = "dR",
  rainbow = FALSE,
  size = 2,
  scale = c(),
  mincol = "yellow",
  maxcol = "red",
  symbol = "feeding",
  symbol.legend = TRUE,
  ocean.col = "aliceblue",
  land.col = rgb(0, 0.5, 0, 0.6)
)
}
\arguments{
\item{longitude}{Longitude of the point. Can only deal with one point at a time.}

\item{latitude}{Latitude of the point. Can only deal with one point at a time.}

\item{nearest}{The number of shell values to be returned. Defaults to 50.}

\item{colour}{The variable to be plotted as colour. Expects a continuous variable. Defaults to 'dR'.}

\item{rainbow}{Whether or not to use a rainbow scale to plot the variable.}

\item{size}{Size of the symbols. Defaults to 2.}

\item{scale}{Resolution of the map. Can be "small", "medium" or "large". If the latter, a high-resolution dataset will have to be downloaded using the R package 'rnaturalearthhires'. Since this package is not on CRAN, you will have to download it yourself. Defaults to 'medium' if 'rnaturalearthhires' is not installed, and to 'high' if it is installed.}

\item{mincol}{Colour for minimum values.}

\item{maxcol}{Colour for maximum values.}

\item{symbol}{The variable to be plotted as symbol. Expects a categoric variable. Defaults to 'feeding'.}

\item{symbol.legend}{Whether or not to plot the legend for the symbols.}

\item{ocean.col}{Colour for the oceans. Defaults to \code{ocean.col="aliceblue"}.}

\item{land.col}{Colour for the land. Defaults to semi-transparent darkgreen: \code{land.col=rgb(0, 0.5, 0, 0.6)}.}
}
\value{
A dataset with the n nearest dR values, and a plot of their coordinates.
}
\description{
Find the shells closest to a chosen coordinate, and plot the dR values and feeding ecology. Uses the marine database downloaded (30 Aug 2024) from calib.org/marine. See Reimer PJ, Reimer RW, 2001. A marine reservoir correction database and on-line interface. Radiocarbon 43:461-3.
}
\examples{
  N_UK <- map.shells(53, -11, 60, 2, scale="medium")
  mean(N_UK$dR)
}
