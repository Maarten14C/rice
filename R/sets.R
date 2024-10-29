# functions to deal with sets of dates

#' @name pool
#' @title Test if a set of radiocarbon dates can be combined 
#' @description Calculate the (chi-square) probability that a set of radiocarbon dates is consistent, i.e. that it can be assumed that they all pertain to the same true radiocarbon age (and thus to the same calendar age - note though that sometimes multiple calendar ages obtain the same C14 age). The function calculates the differences (chi2 value) and finds the corresponding p-value. If the chi2 values is sufficiently small, then the p-value is sufficiently large (above the threshold), and the pooled mean is calculated and returned. If the scatter is too large, no pooled mean is calculated. 
#' @details This follows the calculations of Ward and Wilson (1978; Archaeometry 20: 19-31 <doi:10.1111/j.1475-4754.1978.tb00208.x>) and should only be used for multiple dates that stem from the same sample (e.g., multiple measurements on a single bone). It cannot be used to test if multiple dates from multiple samples pertain to the same event. Since the assumption is that all measurements stem from the same event, we can assume that they all share the same C14 age (since any calBP age will have an associated IntCal C14 age).
#' @return The pooled mean and error if the p-value is above the threshold - a warning if it is not.
#' @param y The set of radiocarbon dates to be tested
#' @param er The lab errors of the radiocarbon dates
#' @param threshold Probability threshold above which chisquare values are considered acceptable (between 0 and 1; default \code{threshold=0.05}).
#' @param roundby Rounding of the reported mean, chisquare and and p-value. Defaults to \code{roundby=1}.
#' @author Maarten Blaauw
#' @examples
#'   data(shroud)
#'   pool(shroud$y,shroud$er)
#'   Zu <- grep("ETH", shroud$ID) # Zurich lab only
#'   pool(shroud$y[Zu],shroud$er[Zu])
#' @export
pool <- function(y, er, threshold=.05, roundby=1) {
  pooled.y <- (sum(y/er^2)) / (sum(1/er^2))
  pooled.er <- sqrt(1/sum(1/er^2))
  T <- sum((y-pooled.y)^2 / er^2)
  p <- pchisq(T, length(y)-1, lower.tail=FALSE)
  
  if(p < threshold) {
    message("! Scatter too large to calculate the pooled mean\nChisq value: ", 
      round(T, roundby+2), ", p-value ", round(p, roundby+2), " < ", threshold, sep="") 
  } else { 
      message("pooled mean: ", round(pooled.y, roundby), " +- ", round(pooled.er, roundby), 
		"\nChi2: ", T, ", p-value ", round(p, roundby+2), " (>", threshold, ", OK)", sep="")
	    return(c(pooled.y, pooled.er))
	  }  
}



#' @name as.one
#' @title Combine multiple radiocarbon dates assuming they belong to the same single year
#' @description Combine all calibrated dates by calculating their product for a range of calendar ages, as if all dates belonged to the same (unknown) single calendar age. This assumed that they all belong to the same single year in time. Use with great care, as often dates could stem from material that could have accumulated over a (much) longer time-span, and if so, then the result will be wrong. See Baillie (1991)'s 'suck-in' effect, Journal of Theoretical Archaeology 2, 12-16. 
#' @details This calculates the product of all calibrated probabilities, over the range of calendar ages to which the radiocarbon ages calibrate. 
#' @return The product of all calibrated probabilities over the range of cal BP years. 
#' @param y The set of radiocarbon dates to be tested
#' @param er The lab errors of the radiocarbon dates
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param is.F Set this to TRUE if the provided age and error are in the F14C realm.
#' @param as.F Whether or not to calculate ages in the F14C realm. Defaults to \code{as.F=FALSE}, which uses the C14 realm.
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
as.one <- function(y, er, cc=1, postbomb=FALSE, is.F=FALSE, as.F=FALSE, thiscurve=NULL, yrsteps=1, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, BCAD=FALSE, cc.dir=NULL, age.lim=c(), age.lab=c(), calib.col=rgb(0,0,0,.2), one.col=rgb(0,0,1,.5), one.height=4, prob=0.95, talk=TRUE, roundby=0, bty="n") {

  calib <- draw.dates(y, er, d.lim=c(0, length(y)+1), BCAD=BCAD, age.lim=age.lim, age.lab=age.lab, col=calib.col, mirror=FALSE, hpd.col=NA, border=NA, yaxt="n", d.lab="", bty=bty)  
  
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

  xpol <- cbind(c(xmin, xseq, xmax), length(y)+1 - c(0, one.height*product/max(product), 0))
  polygon(xpol, col=one.col, border=one.col)
  as.dist <- cbind(xseq, product/max(product))
  as.hpd <- hpd(as.dist, prob, rounded=roundby)
  for(i in 1:nrow(as.hpd))
    segments(as.hpd[i,1], length(y)+1, as.hpd[i,2], length(y)+1, lwd=3, col=rgb(0,0,0,0.5), lend=2)
  
  if(talk) {
    as.points <- suppressWarnings(point.estimates(as.dist, rounded=roundby))
    message("point estimates (mean, median, mode and midpoint): ", as.points[1], ", ", as.points[2], ", ", as.points[3], " & ", as.points[4], ifelse(BCAD, " BC/AD", " cal BP"))

	hpds <- paste0(100 * prob, "% hpd ranges: ", as.hpd[1, 1], "-", as.hpd[1, 2], " (", as.hpd[1, 3], "%)")
	if(nrow(as.hpd) > 1) 
	  for(i in 2:nrow(as.hpd)) 
	    hpds <- paste0(hpds, ", ", as.hpd[i, 1], "-", as.hpd[i, 2], " (", as.hpd[i, 3], "%)")
	message(hpds)
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
#' @param is.F Set this to TRUE if the provided age and error are in the F14C realm.
#' @param as.F Whether or not to calculate ages in the F14C realm. Defaults to \code{as.F=FALSE}, which uses the C14 realm.
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
#' @param calib.col The colour of the individual calibrated ages. Defaults to semi-transparent grey.
#' @param one.col The colour of the combined
#' @param one.height The height of the combined distribution
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
as.bin <- function(y, er, width=100, move.by=c(), move.res=100, cc=1, postbomb=FALSE, is.F=FALSE, as.F=FALSE, thiscurve=NULL, yrsteps=1, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, BCAD=FALSE, cc.dir=NULL, age.lim=c(), age.lab=c(), calib.col=rgb(0,0,0,.2), one.col=rgb(0,0,1,.5), one.height=1, talk=TRUE, prob=0.95, roundby=0, bty="n") {
	
  calib <- draw.dates(y, er, d.lim=c(0, length(y)+1), BCAD=BCAD, age.lim=age.lim, age.lab=age.lab, col=calib.col, mirror=FALSE, hpd.col=NA, border=NA, yaxt="n", d.lab="", bty=bty)
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
  pol <- cbind(c(min(xseq), xseq, max(xseq)), length(y)+1-c(0, one.height*inbin, 0))
  polygon(pol, col=one.col, border=one.col)
  axis(2, length(y)+1-pretty(inbin), labels=pretty(inbin), col=one.col, col.axis=one.col, cex=.8, mgp=c(1.7, .5, 0))
  maxbin <- xseq[which(inbin == max(inbin))[1]]
  abline(v=maxbin+c(-1*width/2, width/2), col=one.col)  

  as.dist <- cbind(xseq, inbin/max(inbin))
  as.hpd <- hpd(as.dist, prob, rounded=roundby)
  for(i in 1:nrow(as.hpd))
    segments(as.hpd[i,1], length(y)+1, as.hpd[i,2], length(y)+1, lwd=3, col=rgb(0,0,0,0.5), lend=2)
  
  if(talk) {
    as.points <- suppressWarnings(point.estimates(as.dist, rounded=roundby))
	message("Most likely bin: ", round(maxbin, roundby), " (", round(maxbin+(width/2), roundby), "-", round(maxbin-(width/2), roundby), " cal BP), fitting a sum of ", round(max(inbin),1), " dates")   
	message("point estimates (mean, median, mode and midpoint): ", as.points[1], ", ", as.points[2], ", ", as.points[3], " & ", as.points[4], ifelse(BCAD, " BC/AD", " cal BP"))

	hpds <- paste0(100 * prob, "% hpd ranges: ", as.hpd[1, 1], "-", as.hpd[1, 2], " (", as.hpd[1, 3], "%)")
	if(nrow(as.hpd) > 1) 
	  for(i in 2:nrow(as.hpd)) 
	    hpds <- paste0(hpds, ", ", as.hpd[i, 1], "-", as.hpd[i, 2], " (", as.hpd[i, 3], "%)")
	message(hpds)
  }
	
  invisible(cbind(xseq, inbin))
}	



#' @name spread
#' @title The spread among calibrated dates
#' @description Calculates the spread among multiple calibrated radiocarbon dates. It does this by randomly sampling ages from the calibrated dates, and calculate the difference between one random date and all others for that iteration. 
#' @return The spread of all calibrated probabilities. 
#' @param y The set of radiocarbon dates
#' @param er The lab errors of the radiocarbon dates
#' @param n The number of iterations to base the calculations on. Defaults to 100,000.
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param as.F Whether or not to calculate ages in the F14C realm. Defaults to \code{as.F=FALSE}, which uses the C14 realm.
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
#' @author Maarten Blaauw
#' @examples
#'   data(shroud)
#'   spread(shroud$y,shroud$er)
#'   Zu <- grep("ETH", shroud$ID) # Zurich lab only
#'   spread(shroud$y[Zu],shroud$er[Zu])
#' @export
spread <- function(y, er, n=1e5, cc=1, postbomb=FALSE, as.F=FALSE, thiscurve=NULL, yrsteps=1, cc.resample=FALSE, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, cc.dir=NULL, visualise=TRUE, talk=TRUE, prob=0.95, roundby=1, bty="l") {
  xs <- array(NA, dim=c(n, length(y)))
  ns <- sample(1:length(y), n, replace=TRUE)
  diffs <- array(NA, dim=c(n, length(y)-1))
  
  for(i in 1:length(y))
    xs[,i] <- r.calib(n, y[i], er[i], cc=cc, postbomb=postbomb, as.F=as.F, thiscurve=thiscurve, yrsteps=yrsteps, cc.resample=cc.resample, threshold=threshold, normal=normal, t.a=t.a, t.b=t.b, cc.dir=cc.dir)
  for(i in 1:n)
    diffs[i,] <- abs(xs[i,ns[i]] - xs[i,-ns[i]])
  
  as.dens <- density(diffs)
  as.hpd <- hpd(diffs)
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
