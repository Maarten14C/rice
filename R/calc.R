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



# internal plotting function
plot_contamination <- function(trueF, obsF, percentage, contamF, ylim=c(), xlab="contamination (%)", true.col="black", observed.col="blue", contamination.col="red", true.pch=20, observed.pch=18, contamination.pch=17, true.name="target", ylab="F14C", bty="l") {
  if(length(ylim) == 0)
    ylim <- sort(extendrange(c(obsF, trueF, contamF)))
  plot(0, type="n", xlim=c(0, 100), ylim=ylim, xlab=xlab, ylab=ylab, bty=bty)
  segments(-100, trueF, 0, trueF, lty=3, col=true.col) # horizontal to target
  segments(0, trueF, 0, -9999, lty=3, col=true.col) # vertical to target
  segments(-100, obsF, percentage, lty=3, col=observed.col) # horizontal from observed
  segments(percentage, obsF, percentage, -9999, lty=3, col=observed.col) # vertical from observed
  segments(-100, contamF, 100, contamF, lty=3, col=2) # horizontal from contamination
  segments(100, contamF, 100, -9999, lty=3, col=2) # vertical from contamination
  segments(0, trueF, 100, contamF, lty=3, col=grey(0.5)) # the connecting line
  points(0, trueF, col=true.col, pch=true.pch)
  points(percentage, obsF, col=observed.col, pch=observed.pch)
  points(100, contamF, col=contamination.col, pch=contamination.pch)

  legend("right", pch=c(true.pch, observed.pch, contamination.pch), col=c(true.col, observed.col, contamination.col), legend=c(true.name, "observed", "contamination"), bg=rgb(1,1,1,0.8), box.lty=0, cex=.7, xjust=0, yjust=0)
}



#' @name contaminate
#' @title Simulate the impact of contamination on a radiocarbon age
#' @description Given a true/target radiocarbon age, calculate the impact of contamination (for example, 1\% contamination with modern carbon) on the observed age
#' @return The observed radiocarbon age and error
#' @param y the true radiocarbon age
#' @param er the error of the true radiocarbon age
#' @param percentage Relative amount of contamination. Must be between 0 and 1
#' @param F.contam the F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param contam.er error of the contamination. Defaults to 0.
#' @param decimals Rounding of the output. Since details matter here, the default is to provide 5 decimals.
#' @param visualise By default, a plot is made to visualise the target and observed F14C values, together with the inferred contamination.
#' @param talk Whether or not to report the calculations made. Defaults to \code{talk=TRUE}.
#' @param true.col Colour for the target/true values. Defaults to black.
#' @param observed.col Colour for the observed values. Defaults to blue.
#' @param contamination.col Colour for the contamination values. Defaults to red.
#' @param true.pch Icon for the true/target date. Defaults to a filled circle.
#' @param observed.pch Icon for the observed. Defaults to a diamond.
#' @param contamination.pch Icon for the contamination. Defaults to a triangle.
#' @param true.name Name of the label of the true/target date
#' @param xlab Name of the x-axis. Defaults to 'contamination (\%)'.
#' @param ylab Name of the y-axis. Defaults to 'F14C'.
#' @param ylim Limits of the y-axis. Calculated automatically by default.
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="l"}.
#' @author Maarten Blaauw
#' @examples
#' contaminate(5000, 20, 1, 1) # 1% contamination with modern carbon
#' contaminate(66e6, 1e6, 1, 1) # dino bone, shouldn't be dated as way beyond the dating limit
#' @export
contaminate <- function(y, er=0, percentage, F.contam=1, contam.er=0, decimals=5, visualise=TRUE, talk=TRUE, true.col="black", observed.col="blue", contamination.col="red", true.pch=20, observed.pch=18, contamination.pch=17, true.name="true", xlab="contamination (%)", ylab="F14C", ylim=c(), bty="l") {
  fraction <- percentage/100
  true.F <- as.data.frame(C14toF14C(y, er, decimals))
  obs.F <- ((1-fraction)*true.F[,1]) + (fraction*F.contam)
  er <- sqrt(true.F[,2]^2 + contam.er^2)
  obs.C14 <- F14CtoC14(obs.F, er, decimals)

  if(visualise)
    if(length(y) == 1)
      plot_contamination(true.F[,1], obs.F, percentage, F.contam, ylim=ylim, xlab=xlab, true.col=true.col, observed.col=observed.col, contamination.col=contamination.col, true.pch=true.pch, true.name=true.name, observed.pch=observed.pch, contamination.pch=contamination.pch, ylab=ylab, bty=bty)
  
  if(talk)
    if(length(y) == 1) { # only report if one value presented
      message("True age as F14C: ", round(true.F[,1], decimals), " +- ", round(true.F[,2], decimals))
      message("Observed F14C: (", 1-fraction, "*", round(true.F[,1],decimals), ") + (", fraction, "*", round(F.contam, decimals), ") = ", round(obs.F, decimals), " +- ", round(er, decimals))
      message("Observed C14 age: ", obs.C14[1], " +- ", obs.C14[2])
    }
  invisible(cbind(obs.C14)) 
}



#' @name clean
#' @title Simulate removing contamination from a radiocarbon age
#' @description Given an observed radiocarbon age, remove the impact of contamination (for example, 1\% contamination with modern carbon) to estimate the true/target age
#' @return The true/target radiocarbon age and error
#' @param y the observed radiocarbon age
#' @param er the error of the observed radiocarbon age
#' @param percentage Relative amount of contamination. Must be between 0 and 100 (\%)
#' @param F.contam the F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param contam.er error of the contamination. Defaults to 0.
#' @param decimals Rounding of the output. Since details matter here, the default is to provide 5 decimals.
#' @param visualise By default, a plot is made to visualise the target and observed F14C values, together with the inferred contamination.
#' @param talk Whether or not to report the calculations made. Defaults to \code{talk=TRUE}.
#' @param true.col Colour for the true/target values. Defaults to black.
#' @param observed.col Colour for the observed values. Defaults to blue.
#' @param contamination.col Colour for the contamination values. Defaults to red.
#' @param true.pch Icon for the true/target date. Defaults to a filled circle.
#' @param observed.pch Icon for the observed. Defaults to a diamond
#' @param contamination.pch Icon for the contamination. Defaults to a triangle.
#' @param true.name Name of the label of the true/target date
#' @param xlab Name of the x-axis. Defaults to 'contamination (\%)'.
#' @param ylab Name of the y-axis. Defaults to 'F14C'.
#' @param ylim Limits of the y-axis. Calculated automatically by default.
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="l"}.
#' @author Maarten Blaauw
#' @examples
#' clean(5000, 20, 1, 1) # 1% contamination with modern carbon
#' @export
clean <- function(y, er=0, percentage, F.contam=1, contam.er=0, decimals=5, visualise=TRUE, talk=TRUE, true.col="black", observed.col="blue", contamination.col="red", true.pch=20, observed.pch=18, contamination.pch=17, true.name="true", xlab="contamination (%)", ylab="F14C", ylim=c(), bty="l") {
  if(length(y)>1)
    stop("cannot deal with more than one value at a time")
  fraction <- percentage/100
  obs.F <- as.data.frame(C14toF14C(y, er, decimals))
  true.F <- (obs.F[,1] - fraction * F.contam) / (1 - fraction)
  er <- sqrt(obs.F[,2]^2 + contam.er^2)
  true.C14 <- F14CtoC14(true.F, er, decimals)

  if(visualise)
    if(length(y) == 1) # not if multiple entries
      plot_contamination(true.F, obs.F[,1], percentage, F.contam, ylim=ylim, xlab=xlab, true.col=true.col, observed.col=observed.col, contamination.col=contamination.col, true.pch=true.pch, true.name=true.name, observed.pch=observed.pch, contamination.pch=contamination.pch, ylab=ylab, bty=bty)
  
  if(talk)
    if(length(y) == 1) {
      message("Observed age as F14C: ", round(obs.F[1], decimals), " +- ", round(obs.F[2], decimals))
      message("True F14C: (", 1-fraction, "*", round(obs.F[1],decimals), ") - (", fraction, "*", round(F.contam, decimals), ") = ", round(true.F, decimals), " +- ", round(er, decimals))
      message("True C14 age: ", true.C14[1], " +- ", true.C14[2])
    }
  invisible(cbind(true.C14)) 
}



#' @name muck
#' @title Calculate the amount of muck/contamination to explain an observed C14 age
#' @description Given an observed and a target radiocarbon age, calculate the amount of contamination required to explain the observed age.
#' @return The required contamination (as percentage), as well as a plot
#' @param y.obs the observed radiocarbon age
#' @param y.target the target radiocarbon age
#' @param F.contam the F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param decimals Rounding of the output. Since details matter here, the default is to provide 5 decimals.
#' @param talk Whether or not to report the calculations made
#' @param visualise By default, a plot is made to visualise the target and observed F14C values, together with the inferred contamination.
#' @param talk Whether or not to report the calculations made. Defaults to \code{talk=TRUE}.
#' @param target.col Colour for the target values. Defaults to black.
#' @param observed.col Colour for the observed values. Defaults to blue.
#' @param contamination.col Colour for the contamination values. Defaults to red.
#' @param target.pch Icon for the target. Defaults to a filled circle.
#' @param observed.pch Icon for the observed. Defaults to a diamond
#' @param contamination.pch Icon for the contamination. Defaults to a triangle.
#' @param true.name Name of the label of the true/target date
#' @param xlab Name of the x-axis. Defaults to 'contamination (\%)'.
#' @param ylab Name of the y-axis. Defaults to 'F14C'.
#' @param ylim Limits of the y-axis. Calculated automatically by default.
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="l"}.
#' @author Maarten Blaauw
#' @examples
#'   muck(600, 2000, 1)
#' @export
muck <- function(y.obs, y.target, F.contam=1, decimals=3, visualise=TRUE, talk=TRUE, target.col="black", observed.col="blue", contamination.col="red", target.pch=20, observed.pch=18, contamination.pch=17, true.name="target", xlab="contamination (%)", ylab="F14C", ylim=c(), bty="l") {
  obsF <- C14toF14C(y.obs)
  targetF <- C14toF14C(y.target)

  # note: these calculations do NOT take into account uncertainties
  frac <- ((obsF - targetF) / (F.contam - targetF))
  perc <- 100*frac

  if(visualise)
    if(length(y.obs) == 1)
      plot_contamination(targetF, obsF, perc, F.contam, ylim=ylim, xlab=xlab, true.col=target.col, observed.col=observed.col, contamination.col=contamination.col, true.pch=target.pch, true.name=true.name, observed.pch=observed.pch, contamination.pch=contamination.pch, ylab=ylab, bty=bty)
 
  obsF <- round(obsF, decimals); targetF <- round(targetF, decimals); perc <- round(perc,decimals)

  if(talk)
    if(length(y.obs) == 1) {
      message("Observed age: ", y.obs, " C14 BP (", obsF, " F14C)")
      message("Target age: ", y.target, " C14 BP (", targetF, " F14C)")
      message("Calculation: (", obsF, "-", targetF, ")/(", F.contam, " - ",  targetF, ") = ", round(frac, decimals+2))
      message("Contamination required: ", perc, "%" )
    }
  invisible(frac)
}
