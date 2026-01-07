# functions to deal with sets of dates

#' @name pool
#' @title Test if a set of radiocarbon dates can be combined 
#' @description Calculate the (chi-square) probability that a set of radiocarbon dates is consistent, i.e. that it can be assumed that they all pertain to the same true radiocarbon age (and thus to the same calendar age - note though that sometimes multiple calendar ages obtain the same C14 age). The function calculates the differences (chi2 value) and finds the corresponding p-value. If the chi2 values is sufficiently small, then the p-value is sufficiently large (above the threshold), and the pooled mean is calculated and returned. If the scatter is too large, no pooled mean is calculated. 
#' @details This follows the calculations of Ward and Wilson (1978; Archaeometry 20: 19-31 <doi:10.1111/j.1475-4754.1978.tb00208.x>) and should only be used for multiple dates that stem from the same sample (e.g., multiple measurements on a single bone). It cannot be used to test if multiple dates from multiple samples pertain to the same event. Since the assumption is that all measurements stem from the same event, we can assume that they all share the same C14 age (since any calBP age will have an associated IntCal C14 age).
#' @return The pooled mean and error if the p-value is above the threshold - a warning if it is not.
#' @param y The set of radiocarbon dates to be tested
#' @param er The lab errors of the radiocarbon dates
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param threshold Probability threshold above which chisquare values are considered acceptable (between 0 and 1; default \code{threshold=0.05}).
#' @param roundby Rounding of the reported mean, chisquare and and p-value. Defaults to \code{roundby=1}.
#' @param talk It's better than staying silent.
#' @author Maarten Blaauw
#' @examples
#'   data(shroud)
#'   pool(shroud$y,shroud$er)
#'   Zu <- grep("ETH", shroud$ID) # Zurich lab only
#'   pool(shroud$y[Zu],shroud$er[Zu])
#' @export
pool <- function(y, er, deltaR=0, deltaSTD=0, threshold=.05, roundby=1, talk=TRUE) {

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  pooled.y <- (sum(y/er^2)) / (sum(1/er^2))
  pooled.er <- sqrt(1/sum(1/er^2))
  T <- sum((y-pooled.y)^2 / er^2)
  p <- pchisq(T, length(y)-1, lower.tail=FALSE)
  
  if(p < threshold) {
    if(talk)
      message("! Scatter too large to calculate the pooled mean\nChisq value: ",
        round(T, roundby+2), ", p-value ", round(p, roundby+2), " < ", threshold, sep="")
  } else { 
      if(talk)
        message("pooled mean: ", round(pooled.y, roundby), " +- ", round(pooled.er, roundby),
         "\nChi2: ", round(T, roundby+2), ", p-value ", round(p, roundby+2), " (>", threshold, ", OK)", sep="")
      return(c(pooled.y, pooled.er))
    }
}



#' @name as.one
#' @title Combine multiple radiocarbon dates assuming they belong to the same single year
#' @description Combine all calibrated dates by calculating their product for a range of calendar ages, as if all dates belonged to the same (unknown) single calendar age. This assumes that they all belong to the same single year in time. Use with great care, as often dates could stem from material that could have accumulated over a (much) longer time-span, and if so, then the result will be wrong. See Baillie (1991)'s 'suck-in' effect, Journal of Theoretical Archaeology 2, 12-16. 
#' @details This calculates the product of all calibrated probabilities, over the range of calendar ages to which the radiocarbon ages calibrate. 
#' @return The product of all calibrated probabilities over the range of cal BP years. 
#' @param y The set of radiocarbon dates to be tested
#' @param er The lab errors of the radiocarbon dates
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param is.F Set this to TRUE if the provided age and error are in the F14C timescale.
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param yrsteps Steps to use for interpolation. Defaults to the cal BP steps in the calibration curve
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=1e-6}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param age.lim Limits of the age axis. Calculated automatically by default.
#' @param age.lab Label of the age axis. Defaults to cal BP or cal BC/AD.
#' @param d.lim Limits of the depth/vertical axis. Calculated automatically by default.
#' @param calib.col The colour of the individual calibrated ages. Defaults to semi-transparent grey.
#' @param one.col The colour of the combined
#' @param one.height The height of the combined distribution
#' @param prob Probability range for highest posterior density (hpd) values. Defaults to \code{prob=0.95}.
#' @param talk Whether or not to provide an analysis of the results
#' @param roundby Rounding of reported years. Defaults to 0 decimals
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="n"}.
#' @author Maarten Blaauw
#' @examples
#'   data(shroud)
#'   as.one(shroud$y,shroud$er, BCAD=TRUE) # but note the scatter!
#'   Zu <- grep("ETH", shroud$ID) # Zurich lab only
#'   as.one(shroud$y[Zu],shroud$er[Zu], BCAD=TRUE)
#' @export
as.one <- function(y, er, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, is.F=FALSE, as.F=FALSE, thiscurve=NULL, yrsteps=1, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, BCAD=FALSE, cc.dir=NULL, age.lim=c(), age.lab=c(), d.lim=c(), calib.col=rgb(0,0,0,.2), one.col=rgb(0,0,1,.5), one.height=.3, prob=0.95, talk=TRUE, roundby=0, bty="n") {

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  if(length(d.lim) == 0)
    d.lim <- c(0, length(y)+1)   
  calib <- draw.dates(y, er, d.lim=d.lim, BCAD=BCAD, age.lim=age.lim, age.lab=age.lab, col=calib.col, mirror=FALSE, up=TRUE, hpd.col=calib.col, border=NA, yaxt="n", d.lab="", bty=bty, prob=prob)  
  
  xmin <- c(); xmax <- c()
  for(i in 1:length(y))	{
    yr <- calib$ages[,i] 
    xmin <- min(xmin, yr)
    xmax <- max(xmax, yr)
  }

  xseq <- seq(xmin, xmax, by=yrsteps)
  probs <- array(NA, dim=c(length(y), length(xseq)))
  for(i in 1:length(y))
    probs[i,] <- approx(calib$ages[,i], calib$probs[,i], xseq, rule=2)$y 
  product <- apply(probs, 2, prod)

  as.dist <- cbind(xseq, product/max(product))
  hpds <- draw.dist(as.dist, y.pos=max(d.lim), prob=prob, dist.col=one.col, fraction=one.height, as.unit=FALSE)
  if(talk) {
    as.points <- suppressWarnings(point.estimates(as.dist, rounded=roundby))
    message("point estimates (mean, median, mode and midpoint): ", as.points[1], ", ", as.points[2], ", ", as.points[3], " & ", as.points[4], ifelse(BCAD, " cal BC/AD", " cal BP"))
    myhpds <- paste0(100 * prob, "% hpd ranges: ", hpds[1, 1], "-", hpds[1, 2], " (", hpds[1, 3], "%)")
    if(nrow(hpds) > 1)
      for(i in 2:nrow(hpds))
        myhpds <- paste0(myhpds, ", ", hpds[i, 1], "-", hpds[i, 2], " (", hpds[i, 3], "%)")
    message(myhpds)
  }
  
  invisible(as.dist)
}



#' @name as.bin
#' @title Combine multiple radiocarbon dates within bins
#' @description Combine all calibrated dates by calculating their product for a range of calendar ages, as if all dates belonged to the same (unknown) calendar age bin. 
#' @details This calculates the amount of calibrated dates that fall within a specific bin, and calculates these bins as moving windows over the range of calendar ages to which the radiocarbon ages calibrate. 
#' @return The number of dates that fall within the moving bins, for each bin.
#' @param y The set of radiocarbon dates to be tested
#' @param er The lab errors of the radiocarbon dates
#' @param width The bin width to apply. Narrower bins will result in fewer dates fitting those bins, but in more detailed bin width histograms. 
#' @param move.by Step size by which the window moves. Left empty by default, and then the moves are set by the parameter move.res.
#' @param move.res The amount of steps taken to make the histogram. Defaults to \code{move.res=100} - a compromise between detail obtained and calculation speed. 
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param is.F Set this to TRUE if the provided age and error are in the F14C timescale.
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param yrsteps Steps to use for interpolation. Defaults to the cal BP steps in the calibration curve
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=1e-6}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param age.lim Limits of the age axis. Calculated automatically by default.
#' @param age.lab Label of the age axis. Defaults to cal BP or BC/AD.
#' @param d.lim Limits of the depth/vertical axis. Calculated automatically by default.
#' @param calib.col The colour of the individual calibrated ages. Defaults to semi-transparent grey.
#' @param bin.col The colour of the combined
#' @param bin.height The height of the combined distribution
#' @param talk Whether or not to report the calculations made. Defaults to \code{talk=TRUE}.
#' @param prob Probability range for highest posterior density (hpd) values. Defaults to \code{prob=0.95}.
#' @param roundby Rounding of reported years. Defaults to 0 decimals
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="n"}.
#' @author Maarten Blaauw
#' @examples
#' \donttest{
#'   data(shroud)
#'   shroudbin <- as.bin(shroud$y, shroud$er, 50, 10) 
#'   # bins of 50 yr, moving by 10 yr, slow
#' }
#' @export
as.bin <- function(y, er, width=100, move.by=c(), move.res=100, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, is.F=FALSE, as.F=FALSE, thiscurve=NULL, yrsteps=1, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, BCAD=FALSE, cc.dir=NULL, age.lim=c(), age.lab=c(), d.lim=c(), calib.col=rgb(0,0,0,.2), bin.col=rgb(0,0,1,.5), bin.height=.3, talk=TRUE, prob=0.95, roundby=0, bty="n") {

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  if(length(d.lim) == 0)
    d.lim <- c(0, length(y)+1)   
  calib <- draw.dates(y, er, d.lim=d.lim, BCAD=BCAD, age.lim=age.lim, age.lab=age.lab, col=calib.col, mirror=FALSE, up=TRUE, hpd.col=calib.col, border=NA, yaxt="n", d.lab="", bty=bty, prob=prob)  
  
  if(length(move.by) == 0)
    xseq <- seq(min(calib$ages)-width, max(calib$ages)+width, length=move.res) else
      xseq <- seq(min(calib$ages)-width, max(calib$ages)+width, by=move.by) 
  xmin <- xseq-(width/2); xmax <- xseq+(width/2)
  
  tmp <- sapply(1:length(y), function(i) {
    sapply(1:length(xseq), function(j) {
      p.range(xmin[j], xmax[j], y[i], er[i], cc=cc, postbomb=postbomb, normal=normal, as.F=as.F, t.a=t.a, t.b=t.b, BCAD=BCAD, threshold=threshold)
    })
  })

  inbin <- rowSums(tmp)
  as.dist <- cbind(xseq, inbin/max(inbin))
  hpds <- draw.dist(as.dist, y.pos=max(d.lim), prob=prob, dist.col=bin.col, as.unit=FALSE, fraction=bin.height)
  
  if(talk) {
    as.points <- suppressWarnings(point.estimates(as.dist, rounded=roundby))
    message("point estimates (mean, median, mode and midpoint): ", as.points[1], ", ", as.points[2], ", ", as.points[3], " & ", as.points[4], ifelse(BCAD, " cal BC/AD", " cal BP"))
    myhpds <- paste0(100 * prob, "% hpd ranges: ", hpds[1, 1], "-", hpds[1, 2], " (", hpds[1, 3], "%)")
    if(nrow(hpds) > 1)
      for(i in 2:nrow(hpds))
        myhpds <- paste0(myhpds, ", ", hpds[i, 1], "-", hpds[i, 2], " (", hpds[i, 3], "%)")
    message(myhpds)
  }

  invisible(as.dist)
}


#' @name hpd.overlap
#' @title Check whether hpds of two distributions overlap
#' @description Checks whether any of the highest posterior densities (hpds) of two distributions overlap. 
#' @return TRUE if at least one of the hpds of distA overlaps with that of distB. 
#' @param distA Distribution A. Expects two columns: values and their probabilities (e.g., caldist(130,10, cc=1)).
#' @param distB Distribution B. Expects two columns: values and their probabilities (e.g., caldist(130,10, cc=1)).
#' @param prob The probability of the highest posterior densities. Defaults to 95\%.
#' @examples
#'   distA <- caldist(130, 20, cc=0) # normal distribution
#'   distB <- caldist(130, 20, cc=1, bombalert=FALSE) # calibrated distribution
#'   plot(distB, type="l")
#'   lines(distA, col=2)
#'   hpd.overlap(distA, distB)
#' @export
hpd.overlap <- function(distA, distB, prob=.95) {
  hpdA <- hpd(distA, prob=prob)	
  hpdB <- hpd(distB, prob=prob)	

  for(i in 1:nrow(hpdA)) 
    for(j in 1:nrow(hpdB)) 
      if(max(hpdA[i,1:2]) >= min(hpdB[j,1:2]) && 
        min(hpdA[i,1:2]) <= max(hpdB[j,1:2]))
          return(TRUE)

  return(FALSE) 
}



#' @name overlap
#' @title The overlap between calibrated C14 dates
#' @description Calculates the amount of overlap (as percentage) between two or more calibrated radiocarbon dates. It does this by taking a sequence of calendar dates 'x' and for each calendar date find the calibrated distribution with the minimum height - this minimum height is taken as the overlap between the dates for that age. This is repeated for all 'x'. The sum of these heights is the overlap, which can reach values from 0 to 100\%. 
#' @return The overlap between all calibrated probabilities as percentage, and a plot. 
#' @param y The set of radiocarbon dates. Alternatively, existing distributions can be provided as a list of distributions, e.g. already-calibrated distributions or distributions derived from age-model estimates.
#' @param er The lab errors of the radiocarbon dates
#' @param labels Labels to be printed for the distributions (optional).
#' @param is.F Set this to TRUE if the provided age and error are in the F14C timescale.
#' @param res The resolution to base the calculations on. Defaults to 1000 steps between the minimum and maximum cal BP (these are calculated from the total calendar age range of all calibrated distributions).
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=1e-6}.
#' @param xlim Age limits of the x-axis. Calculated automatically by default.  
#' @param cal.rev Reverse the calendar axis. Defaults to TRUE.
#' @param xlab Label of the calendar age, defaults to BCAD or cal BP.
#' @param yrby Resolution in years. Defaults to \code{by=1}.
#' @param dist.col The colour of the individual (calibrated) distributions. Defaults to semi-transparent grey. Different colours can also be provided for the individual distributions.
#' @param overlap.col The colour of the overlap distribution.
#' @param overlap.border The colour of the border of the overlap distribution.
#' @param overlap.height The height of the overlap distribution.
#' @param talk Whether or not to report a summary of the spread.
#' @param visualise Whether or not to plot the individual distributions and the overlap.
#' @param prob Probability range to report. Defaults to \code{prob=0.95}.
#' @param roundby Number of decimals to report.
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="n"}.
#' @param yaxt Type of y-axis. Defaults to none drawn (\code{yaxt="n"}).
#' @examples
#'   y <- c(3820, 4430) # the C14 ages of a twig and a marine shell from a single layer
#'   er <- c(40, 40) # their lab errors
#'   overlap(y, er, cc=1:2, dist.col=3:4, labels=c("twig", "shell"))
#'   mydists <- list(caldist(130,20, cc=1, bombalert=FALSE), caldist(150, 20, cc=0))
#'   overlap(mydists)
#' @export
overlap <- function(y, er=c(), labels=c(), is.F=FALSE, res=1e3, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, thiscurve=NULL, BCAD=FALSE, normal=TRUE, t.a=3, t.b=4, cc.dir=NULL, threshold=0, xlim=c(), cal.rev=TRUE, xlab=c(), yrby=1, dist.col=rgb(0,0,0,.2), overlap.col=rgb(0,0,1,.4), overlap.border=NA, overlap.height=1, talk=TRUE, visualise=TRUE, prob=0.95, roundby=1, bty="n", yaxt="n") {
  if(length(cc) == 1)
    cc <- rep(cc,length(y))
  if(length(deltaR) == 1)
    deltaR <- rep(deltaR, length(y))
  if(length(deltaSTD) == 1)
    deltaSTD <- rep(deltaSTD, length(y))

  if(!is.list(y)) { # then put y&er into a list
    if(length(y) != length(er))
      stop("length of 'y' should be the same as length of 'er'")

    dists <- list()
    for(i in 1:length(y))
      dists[[i]] <- caldist(y[i], er[i], cc=cc[i], is.F=is.F, postbomb=postbomb,
        deltaR=deltaR[i], deltaSTD=deltaSTD[i], normal=normal, t.a=t.a, t.b=t.b, thiscurve=thiscurve, BCAD=BCAD)
  } else
    dists <- y

  xmin <- Inf; xmax <- -Inf
  for(i in 1:length(dists)) {
    xmin <- min(xmin, dists[[i]][,1])
    xmax <- max(xmax, dists[[i]][,1])
  }

  xseq <- seq(xmin, xmax, by=yrby)
  probs <- array(NA, dim=c(length(xseq), length(dists)))
  for(i in 1:length(dists)) { # find and normalise all distributions
    this.prob <- approx(dists[[i]][,1], dists[[i]][,2], xseq, rule=1)$y
    this.prob[is.na(this.prob)] <- 0
    this.prob <- this.prob / sum(this.prob)
    this.prob[this.prob < threshold] <- 0
    probs[,i] <- this.prob
  }  
  this.overlap <- apply(probs, 1, min, na.rm=TRUE) # find the minimum prob for each x in xseq
  if(length(this.overlap) == 0)
    this.overlap <- 0

  if(length(dist.col) == 1)
    dist.col <- rep(dist.col, length(dists))
  if(length(labels) > 0)
    if(length(labels) < length(dists))
      message("please provide a label for all entries")

  if(visualise) {
    if(length(xlim) == 0)  
      if(BCAD)
        xrng <- range(xseq) else
          xrng <- rev(range(xseq))
    if(!cal.rev)
      xrng <- rev(xrng)	  
    if(length(xlab) == 0)
      if(BCAD)
        xlab <- "cal BC/AD" else
          xlab <- "cal BP"
    plot(0, type="n", xlim=xrng, xlab=xlab, ylim=c(0, length(dists)+1), ylab="", bty=bty, yaxt=yaxt)

    for(i in 1:length(dists)) {
      thisdist <- cbind(dists[[i]][,1], .8*dists[[i]][,2]/max(dists[[i]][,2]))
      draw.dist(thisdist, dist.col=dist.col[i], y.pos=i, dist.border=dist.col[i], prob=prob, mirror=FALSE, BCAD=BCAD)
      if(length(labels) > 0)
        text(max(dists[[i]][,1]), i, labels=labels[i], col=dist.col[i], adj=c(0,-1))
    }
    if(max(this.overlap) > 0)
      draw.dist(cbind(xseq, .8*this.overlap/max(this.overlap)), dist.col=overlap.col, dist.border=overlap.col, prob=prob, y.pos=0, BCAD=BCAD, hpd=FALSE)
  }
  
  if(talk)
    message("overlap: ", round(100*sum(this.overlap)*yrby, roundby), "%")
  invisible(100*sum(this.overlap)*yrby)
}



#' @name spread
#' @title The spread among calibrated dates
#' @description Calculates the spread among multiple calibrated radiocarbon dates. It does this by randomly sampling ages from the calibrated dates, and calculating the difference between one random date and all others for that iteration.
#' @return The spread of all calibrated probabilities. 
#' @param y The set of radiocarbon dates
#' @param er The lab errors of the radiocarbon dates
#' @param n The number of iterations to base the calculations on. Defaults to 100,000. Different values for n could significantly alter performance and accuracy.
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param yrsteps Steps to use for interpolation. Defaults to the cal BP steps in the calibration curve
#' @param cc.resample The IntCal20 curves have different densities (every year between 0 and 5 kcal BP, then every 5 yr up to 15 kcal BP, then every 10 yr up to 25 kcal BP, and then every 20 yr up to 55 kcal BP). If calibrated ages span these density ranges, their drawn heights can differ, as can their total areas (which should ideally all sum to the same size). To account for this, resample to a constant time-span, using, e.g., cc.resample=5 for 5-yr timespans.
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=1e-6}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param visualise Whether or not to plot the spread
#' @param talk Whether or not to report a summary of the spread 
#' @param prob Probability range to report. Defaults to \code{prob=0.95}.
#' @param roundby Number of decimals to report
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="l"}.
#' @examples
#'   data(shroud)
#'   spread(shroud$y,shroud$er)
#'   Zu <- grep("ETH", shroud$ID) # Zurich lab only
#'   spread(shroud$y[Zu],shroud$er[Zu])
#' @export
spread <- function(y, er, n=1e5, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, as.F=FALSE, thiscurve=NULL, yrsteps=1, cc.resample=FALSE, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, cc.dir=NULL, visualise=TRUE, talk=TRUE, prob=0.95, roundby=1, bty="l") {

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  xs <- array(NA, dim=c(n, length(y)))
  ns <- sample(1:length(y), n, replace=TRUE)
  diffs <- array(NA, dim=c(n, length(y)-1))
  
  for(i in 1:length(y))
    xs[,i] <- r.calib(n, y[i], er[i], cc=cc, postbomb=postbomb, as.F=as.F, thiscurve=thiscurve, yrsteps=yrsteps, cc.resample=cc.resample, threshold=threshold, normal=normal, t.a=t.a, t.b=t.b, cc.dir=cc.dir)
  for(i in 1:n)
    diffs[i,] <- abs(xs[i,ns[i]] - xs[i,-ns[i]])
  as.dens <- density(diffs)
  as.hpd <- hpd(diffs, every=min(diff(as.dens$x)))
  aspoints <- round(point.estimates(cbind(as.dens$x, as.dens$y)), roundby)
  oneprob <- (1-prob)/2
  minmax <- round(quantile(diffs, probs=c(oneprob, .5, 1-oneprob)), roundby)
  if(talk) 
    message("average spread: ", aspoints[1], " calendar years (",
      "median ", aspoints[2], ")\n95% range: ", minmax[1], " to ", minmax[3])
  
  if(visualise) {
    plot(density(diffs, from=0), xlab="spread (calendar years)", main="", bty=bty)
    segments(minmax[1], 0, minmax[3], 0, lwd=4, col=rgb(0,0,0,.5), lend=1)
    abline(v=c(mean(diffs), minmax[2]), col=c(2,4), lty=2)
    legend("topright", legend=c("spread", "range", "median", "mean"),
      text.col=c(1,grey(.5), 4, 2), bty="n", cex=.7)
  }
  
  invisible(diffs)
}



#' @name span
#' @title The time span between two calibrated dates
#' @description Calculates the timespan between two calibrated radiocarbon dates. It does this by randomly sampling ages from both calibrated dates, followed by calculating the differences between all samples ages.
#' @return The time span.
#' @param y1 The first radiocarbon date.
#' @param er1 The lab error of the first radiocarbon date.
#' @param y2 The second radiocarbon date.
#' @param er2 The lab error of the second radiocarbon date.
#' @param n The number of iterations to base the calculations on. Defaults to 100,000. Different values for n could significantly alter performance and accuracy.
#' @param positive Whether or not to enforce the span to be positive. If set to TRUE, then negative span values are removed. Defaults to TRUE.
#' @param cc Calibration curve(s) to use. Defaults to IntCal20 (\code{cc=1}). Can be a vector of length 2.
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param deltaR Age offset (e.g. for marine samples). Can be a vector of length 2.
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation). Can be a vector of length 2.
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param yrsteps Steps to use for interpolation. Defaults to the cal BP steps in the calibration curve
#' @param cc.resample The IntCal20 curves have different densities (every year between 0 and 5 kcal BP, then every 5 yr up to 15 kcal BP, then every 10 yr up to 25 kcal BP, and then every 20 yr up to 55 kcal BP). If calibrated ages span these density ranges, their drawn heights can differ, as can their total areas (which should ideally all sum to the same size). To account for this, resample to a constant time-span, using, e.g., cc.resample=5 for 5-yr timespans.
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=1e-6}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3). Can be a vector of length 2.
#' @param t.b Value b of the t distribution (defaults to 4). Can be a vector of length 2.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param visualise Whether or not to plot the time span.
#' @param talk Whether or not to report a summary of the span.
#' @param prob Probability range to report. Defaults to \code{prob=0.95}.
#' @param roundby Number of decimals to report
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="l"}.
#' @examples
#'   span(2300, 20, 2350, 20)
#' @export
span <- function(y1, er1, y2, er2, n=1e5, positive=TRUE, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, as.F=FALSE, thiscurve=NULL, yrsteps=1, cc.resample=FALSE, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, cc.dir=NULL, visualise=TRUE, talk=TRUE, prob=0.95, roundby=1, bty="l") {

  if(!length(y1) == 1)
    stop("y1 must be of length 1")
  if(!length(er1) == 1)
    stop("er1 must be of length 1")
  if(!length(y2) == 1)
    stop("y2 must be of length 1")
  if(!length(er2) == 1)
    stop("er2 must be of length 1")
  if(length(cc) == 1)
    cc <- c(cc, cc)
  if(length(t.a) == 1)
    t.a <- c(t.a, t.a)
  if(length(t.b) == 1)
    t.b <- c(t.b, t.b)
  if(length(deltaSTD) == 1)
    deltaSTD <- c(deltaSTD, deltaSTD)
  if(length(deltaR) == 1)
    deltaR <- c(deltaR, deltaR)
  if(n < 100)
    message("very few iterations, results may be unreliable")
  if(n > 1e6)
    message("very many iterations, overkill, this may take a while and be difficult to compute")

  y1 <- y1 - deltaR[1]
  er1 <- sqrt(er1^2 + deltaSTD[1]^2)
  y2 <- y2 - deltaR[2]
  er2 <- sqrt(er2^2 + deltaSTD[2]^2)

  x1 <- r.calib(n, y1, er1, cc=cc[1], postbomb=postbomb, as.F=as.F, thiscurve=thiscurve, yrsteps=yrsteps, cc.resample=cc.resample, threshold=threshold, normal=normal, t.a=t.a[1], t.b=t.b[1], cc.dir=cc.dir)
  x2 <- r.calib(n, y2, er2, cc=cc[2], postbomb=postbomb, as.F=as.F, thiscurve=thiscurve, yrsteps=yrsteps, cc.resample=cc.resample, threshold=threshold, normal=normal, t.a=t.a[2], t.b=t.b[2], cc.dir=cc.dir)
  diffs <- x2 - x1
  if(positive) {
    diffs <- diffs[which(diffs >= 0)]
    as.dens <- density(diffs, from=0)
    } else
      as.dens <- density(diffs)

  as.hpd <- hpd(cbind(as.dens$x, as.dens$y))
  aspoints <- round(point.estimates(cbind(as.dens$x, as.dens$y)), roundby)
  oneprob <- (1-prob)/2
  minmax <- round(quantile(diffs, probs=c(oneprob, .5, 1-oneprob)), roundby)
  if(talk)
    message("average span: ", aspoints[1], " calendar years (",
      "median ", aspoints[2], ")\n95% range: ", minmax[1], " to ", minmax[3])

  if(visualise) {
    plot(as.dens, xlab="span (calendar years)", main="", bty=bty)
    segments(minmax[1], 0, minmax[3], 0, lwd=4, col=rgb(0,0,0,.5), lend=1)
    abline(v=c(mean(diffs), minmax[2]), col=c(2,4), lty=2)
    legend("topright", legend=c("span", "range", "median", "mean"),
      text.col=c(1,grey(.5), 4, 2), bty="n", cex=.7)
  }

  invisible(diffs)
}
