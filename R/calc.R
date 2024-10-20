#' @name fractions
#' @title Estimate a missing radiocarbon age from fractions
#' @description Estimate a missing radiocarbon age from a sample which has C14 dates on both the bulk and on fractions, but where 1 sample was too small to be dated. This can be used in for example soils separated into size fractions, where one of the samples turns out to be too small to be dated. Requires to have the bulk age, the ages of the dated fractions, and the carbon contents and weights of all fractions.
#' @param bulk_age The age of the bulk/entire sample
#' @param bulk_er The error of the age of the bulk/entire sample
#' @param fractions_percC The \%carbon contents of the fractions. If unknown, enter estimates (e.g., rep(1,4))
#' @param fractions_weights The weights of the fractions. The units are not important here as the weights are used to calculate the relative contributions of carbon within individual fractions to the entire sample.
#' @param fractions_ages The radiocarbon ages of the individual fractions. The fraction without a date should be entered as NA.
#' @param fractions_errors The errors of the radiocarbon ages of the individual fractions. The fraction without a date should be entered as NA.
#' @param roundby Rounding of the reported age
#' @examples
#' Cs <- c(.02, .05, .03, .04) # carbon contents of each fraction
#' wghts <- c(5, 4, 2, .5) # weights for all fractions, e.g., in mg
#' ages <- c(130, 130, 130, NA) # ages of all fractions. The unmeasured one is NA
#' errors <- c(10, 12, 10, NA) # errors, unmeasured is NA
#' fractions(150, 20, Cs, wghts, ages, errors) # assuming a bulk age of 150 +- 20 C14 BP
#' @export
fractions <- function(bulk_age, bulk_er, fractions_percC, fractions_weights, fractions_ages, fractions_errors, roundby=1) {

  if(length(which(is.na(fractions_ages))) > 1 || length(which(is.na(fractions_errors))) > 1)
	stop("Cannot deal with multiple missing fraction ages/errors")  

  unknown_age <- which(is.na(fractions_ages))
  totC <- fractions_percC * fractions_weights # how much C in total
  totC <- totC / sum(totC) # normalise to 1

  # Bulk F14C value and its error
  bulk_F <- C14toF14C(bulk_age, bulk_er)

  # F14C values for the known fractions (ignoring the unknown one)
  fractions_F <- C14toF14C(fractions_ages[-unknown_age], fractions_errors[-unknown_age])

  # Carbon contribution * F14C value for each known fraction
  fractions_cF <- totC[-unknown_age] * fractions_F[,1] # carbon contribution * C14 ages

  # Propagate the uncertainty for the known fractions, weighted by their %C and errors
  total_known_error <- sum((totC[-unknown_age]^2) / (fractions_F[,2]^2))
  overall_uncertainty <- sqrt(total_known_error + bulk_F[2]^2)

  # Calculate the remaining fraction's F14C value and age
  unknown_F <- bulk_F[1] - sum(fractions_cF)
  unknown_age_estimated <- F14CtoC14(unknown_F / totC[unknown_age])

  unknown_age <- round(c(unknown_age_estimated[1], overall_uncertainty), roundby)

  message(paste0("estimated C14 age of fraction ", which(is.na(fractions_ages)), ": ", unknown_age[1], " +- ", unknown_age[2]))
  invisible(unknown_age)
}



#' @name contaminate
#' @title Simulate the impact of contamination on a radiocarbon age
#' @description Given a certain radiocarbon age, calculate the observed impact of contamination with a ratio of material with a different 14C content (for example, 1% contamination with modern carbon)
#' @return The observed radiocarbon age and error
#' @param y the true radiocarbon age
#' @param sdev the error of the true radiocarbon age
#' @param fraction Relative amount of contamination. Must be between 0 and 1
#' @param F14C the F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param F14C.er error of the contamination. Defaults to 0.
#' @param decimals Rounding of the output. Since details matter here, the default is to provide 5 decimals.
#' @author Maarten Blaauw
#' @examples
#' contaminate(5000, 20, .01, 1) # 1% contamination with modern carbon
#' # Impacts of different amounts of contamination with modern carbon:
#' real.14C <- seq(0, 50e3, length=200)
#' contam <- seq(0, .1, length=101) # 0 to 10% contamination
#' contam.col <- rainbow(length(contam))
#' plot(0, type="n", xlim=c(0, 55e3), 
#'   xlab="real", ylim=range(real.14C), ylab="observed")
#' for(i in 1:length(contam))
#'   lines(real.14C, contaminate(real.14C, c(), contam[i], 1, decimals=5), col=contam.col[i])
#' contam.legend <- seq(0, .1, length=6)
#' contam.col <- rainbow(length(contam.legend))
#' text(52e3, contaminate(50e3, c(), contam.legend, 1), labels=contam.legend, col=contam.col, cex=.7)
#' @export
contaminate <- function(y, sdev=c(), fraction, F14C, F14C.er=0, decimals=5) {
  y.F <- as.data.frame(C14toF14C(y, sdev, decimals))
  mn <- ((1-fraction)*y.F[,1]) + (fraction*F14C)
  if(length(sdev) == 0)
    return(F14CtoC14(mn, c(), decimals)) else {
      er <- sqrt(y.F[,2]^2 + F14C.er^2)
      return(F14CtoC14(mn, er, decimals))
    }
}



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
#'   y <- c(280, 346, 359, 387)
#'   y2 <- c(100, 346, 359, 387)
#'   er <- c(15, 13, 16, 18)
#'   pool(y,er)
#'   pool(y2,er)
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


