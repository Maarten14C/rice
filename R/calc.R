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
plot_contamination <- function(true.F, true.er, obs.F, obs.er, perc, perc.er, contam.F, contam.er, ylim=c(), xlab="contamination (%)", true.col="black", observed.col="blue", contamination.col="red", true.pch=20, observed.pch=18, contamination.pch=17, true.name="target", ylab="F14C", bty="l") {
  if(length(ylim) == 0)
    ylim <- sort(extendrange(c(obs.F, true.F, contam.F)))
  plot(0, type="n", xlim=c(0, 100), ylim=ylim, xlab=xlab, ylab=ylab, bty=bty)
  segments(-100, true.F, 0, true.F, lty=3, col=true.col) # horizontal to target
  segments(0, true.F, 0, -9999, lty=3, col=true.col) # vertical to target
  segments(-100, obs.F, perc, lty=3, col=observed.col) # horizontal from observed
  segments(perc, obs.F, perc, -9999, lty=3, col=observed.col) # vertical from observed
  segments(-100, contam.F, 100, contam.F, lty=3, col=2) # horizontal from contamination
  segments(100, contam.F, 100, -9999, lty=3, col=2) # vertical from contamination
  segments(0, true.F, 100, contam.F, lty=3, col=grey(0.5)) # the connecting line

  if(obs.er > 0)
    segments(perc, obs.F-obs.er, perc, obs.F+obs.er, col=observed.col, lwd=2)
  if(perc.er > 0)
    segments(perc-perc.er, obs.F, perc+perc.er, obs.F, col=observed.col, lwd=2)
  if(contam.er > 0)
    segments(100, contam.F+contam.er, 100, contam.F-contam.er, col=contamination.col, lwd=2)
  if(true.er > 0)
    segments(0, true.F-true.er, 0, true.F+true.er, col=true.col, lwd=2)

  points(0, true.F, col=true.col, pch=true.pch)
  points(perc, obs.F, col=observed.col, pch=observed.pch)
  points(100, contam.F, col=contamination.col, pch=contamination.pch)

  legend("right", pch=c(true.pch, observed.pch, contamination.pch), col=c(true.col, observed.col, contamination.col), legend=c(true.name, "observed", "contamination"), bg=rgb(1,1,1,0.8), box.lty=0, cex=.7, xjust=0, yjust=0)
}



#' @name contaminate
#' @title Simulate the impact of contamination on a radiocarbon age
#' @description Given a true/target radiocarbon age, calculate the impact of contamination (for example, 1\% contamination with modern carbon) on the observed age. Can optionally include contamination uncertainties, but then Monte Carlo iterations should be used (option MC=TRUE).
#' @return The observed radiocarbon age and error
#' @param y The 'true' radiocarbon age
#' @param er The error of the 'true' radiocarbon age
#' @param percentage Relative amount of contamination. Must be between 0 and 1
#' @param percentage.error Uncertainty of the contamination. Assumed to be normally distributed (which fails close to 0\% or 100\% contamination levels). Defaults to 0\%.
#' @param F.contam the F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param F.contam.er error of the contamination. Defaults to 0.
#' @param MC Whether or not to use Monte Carlo iterations to estimate the values. Defaults to TRUE, because it treats uncertainties better than if set to FALSE.
#' @param its Amount of iterations to use if MC=TRUE. Defaults to 10,000.
#' @param decimals Rounding of the output for F values. Since details matter here, the default is to provide 5 decimals.
#' @param roundby Rounding of the output for C14 ages. Defaults to 1 decimal.
#' @param visualise By default, a plot is made to visualise the target and observed F14C values, together with the inferred contamination.
#' @param talk Whether or not to report the calculations made. Defaults to \code{talk=TRUE}.
#' @param eq.x Leftmost location of the equation. Can be set to values outside of (0,100) to make the equation invisible.
#' @param eq.y Vertical location of the equation. Defaults to the top of the graph.
#' @param eq.size Size of the font of the equation.
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
#' contaminate(5000, 20, 1, 0, 1) # 1% contamination with modern carbon
#' contaminate(66e6, 1e6, 1, 0, 1) # dino bone, shouldn't be dated as way beyond the dating limit
#' @export
contaminate <- function(y, er=0, percentage, percentage.error=0, F.contam=1, F.contam.er=0, MC=TRUE, its=1e4, decimals=5, roundby=1, visualise=TRUE, talk=TRUE, eq.x=0, eq.y=c(), eq.size=0.7, true.col="black", observed.col="blue", contamination.col="red", true.pch=20, observed.pch=18, contamination.pch=17, true.name="true", xlab="contamination (%)", ylab="F14C", ylim=c(), bty="l") {
  if(percentage < 0 || percentage > 100)
    stop("percentage should be between 0 and 100%", call.=FALSE) 
  if(F.contam < 0)
    stop("F.contam cannot be smaller than 0", call.=FALSE) 

  fraction <- percentage/100
  if(length(er) == 1) {
    if(er == 0)
      er <- rep(0, length(y))
  } else
    if(length(y) != length(er))
	  stop("y and er need to be of the same length", call.=FALSE)

  F.true <- as.data.frame(C14toF14C(y, er, decimals))

  if(MC) {
    if(its < 10)
      stop("too few samples for robust analysis. Increase 'its'", call.=FALSE)
    alpha <- (fraction / (percentage.error / 100))^2
    beta <- ((1 - fraction) / (percentage.error / 100))^2
    if(percentage.error == 0) # then waste no time sampling
	  fraction.samples <- fraction else {
	    fraction.samples <- rbeta(its, shape1 = alpha, shape2 = beta) # samples between 0 and 1
        fraction.samples <- fraction.samples * fraction / mean(fraction.samples) # normalise
      }
    if(all(er == 0)) # then waste no time sampling
      F.true.samples <- F.true[,1] else
        F.true.samples <- rnorm(length(fraction.samples), F.true[,1], F.true[,2])
    if(F.contam.er == 0) # then waste no time sampling
      F.contam.samples <- F.contam else
        F.contam.samples <- rnorm(length(fraction.samples), F.contam, F.contam.er)

    F.obs.samples <- (1 - fraction.samples) * F.true.samples + fraction.samples * F.contam.samples # calculate observed F's
    C14.obs.samples <- F14CtoC14(F.obs.samples, 0, decimals)[,1]
    F.obs <- mean(F.obs.samples)
    F.obs.er <- sd(F.obs.samples)
    if(is.na(F.obs.er)) { # if errors are 0, then sd(0) gives NA
      F.obs.er <- 0
	  C14.obs.er <- 0	
    } else
	   C14.obs.er <- sd(C14.obs.samples)
    C14.obs <- cbind(mean(C14.obs.samples), C14.obs.er)

  } else {
      F.obs <- ((1-fraction)*F.true[,1]) + (fraction*F.contam)
      F.obs.er <- sqrt(F.true[,2]^2 + F.contam.er^2)
      C14.obs <- F14CtoC14(F.obs, er, decimals)
    }

  if(visualise)
    if(length(y) == 1) {
      plot_contamination(true.F=F.true[,1], true.er=F.true[,2], obs.F=F.obs, obs.er=F.obs.er, perc=percentage, perc.er=percentage.error, contam.F=F.contam, contam.er=F.contam.er, ylim=ylim, xlab=xlab, true.col=true.col, observed.col=observed.col, contamination.col=contamination.col, true.pch=true.pch, true.name=true.name, observed.pch=observed.pch, contamination.pch=contamination.pch, ylab=ylab, bty=bty)
      txt <- c("F_obs = ", "(", 1-fraction, "*", round(F.true[,1], decimals), ")", " + ", 
        "(", fraction, "*", round(F.contam, decimals), ")", " = ", round(F.obs, decimals))
      colours <- c(observed.col, rep(true.col, 5), grey(0.5), rep(contamination.col, 5), rep(observed.col, 3))
      xpos <- eq.x+(eq.size*cumsum(c(0, strwidth(txt[-length(txt)]))))
      if(length(eq.y) == 0)
        eq.y <- 1.01*max(c(F.obs, F.true[,1], F.contam))
      for(i in seq_along(txt))
        text(x = xpos[i], y = eq.y, labels = txt[i], col = colours[i], adj = c(0, 0), cex=eq.size)  
    }	  

  C14.obs <- round(C14.obs, roundby)

  if(talk)
    if(length(y) == 1) { # only report if one value presented
      message("True age in F14C: ", round(F.true[,1], decimals), " +- ", round(F.true[,2], decimals))
      message("Observed F14C: (", 1-fraction, "*", round(F.true[,1], decimals), ") + (", fraction, "*", round(F.contam, decimals), ") = ", round(F.obs, decimals), " +- ", round(F.obs.er, decimals))
      message("Observed C14 age: ", C14.obs[1], " +- ", C14.obs[2])
    }
  if(MC)
    invisible(list(obs=C14.obs, F14C=c(mean=F.obs, sd=F.obs.er), samples=cbind(F.true.samples, F.contam.samples, fraction.samples, F.obs.samples))) else
      invisible(cbind(C14.obs))
}



#' @name clean
#' @title Simulate removing contamination from a radiocarbon age
#' @description Given an observed radiocarbon age, remove the impact of contamination (for example, 1\% contamination with modern carbon) to estimate the true/target age
#' @return The true/target radiocarbon age and error
#' @param y The observed radiocarbon age
#' @param er The error of the observed radiocarbon age
#' @param percentage Relative amount of contamination. Must be between 0 and 100 (\%)
#' @param percentage.error Uncertainty of the contamination. Assumed to be normally distributed (which fails close to 0\% or 100\% contamination levels). Defaults to 0\%.
#' @param F.contam The F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param F.contam.er The error of the contamination. Defaults to 0.
#' @param MC Whether or not to use Monte Carlo iterations to estimate the values. Defaults to TRUE, because it treats uncertainties better than if set to FALSE.
#' @param its Amount of iterations to use if MC=TRUE. Defaults to 10,000.
#' @param decimals Rounding of the output. Since details matter here, the default is to provide 5 decimals.
#' @param roundby Rounding of the output for C14 ages. Defaults to 1 decimal.
#' @param visualise By default, a plot is made to visualise the target and observed F14C values, together with the inferred contamination.
#' @param talk Whether or not to report the calculations made. Defaults to \code{talk=TRUE}.
#' @param eq.x Leftmost location of the equation. Can be set to values outside of (0,100) to make the equation invisible.
#' @param eq.y Vertical location of the equation. Defaults to the top of the graph.
#' @param eq.size Size of the font of the equation.
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
#' # 1% contamination with modern carbon (no uncertainties in contamination's percentage or F)
#' clean(5000, 20, 1, 0, 1, 0)
#' # now with errors:
#' clean(5000, 20, 1, 0.1, 1, 0.1)
#' @export
clean <- function(y, er=0, percentage, percentage.error=0, F.contam=1, F.contam.er=0, MC=TRUE, its=1e4, roundby=1, decimals=5, visualise=TRUE, talk=TRUE, eq.x=0, eq.y=c(), eq.size=0.7, true.col="black", observed.col="blue", contamination.col="red", true.pch=20, observed.pch=18, contamination.pch=17, true.name="true", xlab="contamination (%)", ylab="F14C", ylim=c(), bty="l") {
  if(length(y)>1)
    stop("cannot deal with more than one value at a time")
  if(percentage < 0 || percentage > 100)
    stop("percentage should be between 0 and 100%", call.=FALSE) 
  if(F.contam < 0)
    stop("F.contam cannot be smaller than 0", call.=FALSE) 

  fraction <- percentage/100
  if(length(er) == 1) {
    if(er == 0)
      er <- rep(0, length(y))
  } else
    if(length(y) != length(er))
      stop("y and er need to be of the same length", call.=FALSE)
  
  F.obs <- as.data.frame(C14toF14C(y, er, decimals))

  if(MC) {
    if(its < 10)
      stop("too few samples for robust analysis. Increase 'its'", call.=FALSE)
    alpha <- (fraction / (percentage.error / 100))^2
    beta <- ((1 - fraction) / (percentage.error / 100))^2
    if(percentage.error == 0) # then waste no time sampling
	  fraction.samples <- fraction else {
	    fraction.samples <- rbeta(its, shape1 = alpha, shape2 = beta) # samples between 0 and 1
        fraction.samples <- fraction.samples * fraction / mean(fraction.samples) # normalise
      }
    if(all(er == 0)) # then waste no time sampling
      F.obs.samples <- F.obs[,1] else	  
	    F.obs.samples <- rnorm(length(fraction.samples), F.obs[,1], F.obs[,2])
    if(F.contam.er == 0) # then waste no time sampling
      F.contam.samples <- F.contam else
        F.contam.samples <- rnorm(length(fraction.samples), F.contam, F.contam.er)

    F.true.samples <- (F.obs.samples - fraction.samples * F.contam.samples) / (1 - fraction.samples)
    C14.true.samples <- F14CtoC14(F.true.samples, 0, decimals)[,1]
    F.true <- mean(F.true.samples)
    F.true.er <- sd(F.true.samples)
    if(is.na(F.true.er)) { # if errors are 0, then sd(0) gives NA
      F.true.er <- 0
	  C14.true.er <- 0	
    } else
	   C14.true.er <- sd(C14.true.samples)
       C14.true <- cbind(mean(C14.true.samples), C14.true.er)
  } else {
      F.true <- (F.obs[,1] - fraction * F.contam) / (1 - fraction)
      F.true.er <- sqrt(F.obs[,2]^2 + F.contam.er^2)
      C14.true <- F14CtoC14(F.true, er, decimals)
    }

  C14.true <- round(C14.true, roundby)

  if(visualise)
    if(length(y) == 1) { # not if multiple entries
      plot_contamination(true.F=F.true, true.er=F.true.er, obs.F=F.obs[,1], obs.er=F.obs[,2], perc=percentage, perc.er=percentage.error, contam.F=F.contam, contam.er=F.contam.er, ylim=ylim, xlab=xlab, true.col=true.col, observed.col=observed.col, contamination.col=contamination.col, true.pch=true.pch, true.name=true.name, observed.pch=observed.pch, contamination.pch=contamination.pch, ylab=ylab, bty=bty)
      txt <- c("F_true = ", "(", 1-fraction, "*", round(F.obs[1], decimals), ")", " - ", 
        "(", fraction, "*", round(F.contam, decimals), ")", " = ", round(F.true, decimals))
      colours <- c(true.col, rep(observed.col, 5), grey(.5), rep(contamination.col, 5), rep(true.col, 3))
      xpos <- eq.x+(eq.size*cumsum(c(0, strwidth(txt[-length(txt)]))))
      if(length(eq.y) == 0)
        eq.y <- 1.01*max(c(F.obs[,1], F.true, F.contam))
      for(i in seq_along(txt))
        text(x = xpos[i], y = eq.y, labels = txt[i], col = colours[i], adj = c(0, 0), cex=eq.size)
    }

  if(talk)
    if(length(y) == 1) {
      message("Observed age as F14C: ", round(F.obs[1], decimals), " +- ", round(F.obs[2], decimals))
      message("True F14C: (", 1-fraction, "*", round(F.obs[1], decimals), ") - (", fraction, "*", round(F.contam, decimals), ") = ", round(F.true, decimals), " +- ", round(F.true.er, decimals))
      message("True C14 age: ", C14.true[1], " +- ", C14.true[2])
    }
  if(MC)
    invisible(list(obs=C14.true, F14C=c(mean=F.true, sd=F.true.er), samples=cbind(F.obs.samples, F.contam.samples, fraction.samples, F.true.samples))) else
      invisible(cbind(C14.true))
}



#' @name muck
#' @title Calculate the amount of muck/contamination to explain an observed C14 age
#' @description Given an observed and a target radiocarbon age, calculate the amount of contamination required to explain the observed age.
#' @return The required contamination (as percentage), as well as a plot
#' @param y.obs The observed radiocarbon age
#' @param y.obs.er The error of the observed radiocarbon age
#' @param y.target the target radiocarbon age
#' @param y.target.er The error of the target radiocarbon age. Not taken into account in the calculations.
#' @param F.contam the F14C of the contamination. Set at 1 for carbon of modern radiocarbon age, at 0 for 14C-free carbon, or anywhere inbetween.
#' @param F.contam.er The error of the contamination. Defaults to 0.
#' @param MC Whether or not to use Monte Carlo iterations to estimate the values. Defaults to TRUE, because it treats uncertainties better than if set to FALSE.
#' @param its Amount of iterations to use if MC=TRUE. Defaults to 10,000.
#' @param filter If set to true (default), iterations that result in fractions <0 or >1 are removed.
#' @param roundby Rounding of the output for C14 ages. Defaults to 1 decimal.
#' @param decimals Rounding of the output. Since details matter here, the default is to provide 5 decimals.
#' @param visualise By default, a plot is made to visualise the target and observed F14C values, together with the inferred contamination.
#' @param talk Whether or not to report the calculations made
#' @param eq.x Leftmost location of the equation. Can be set to values outside of (0,100) to make the equation invisible.
#' @param eq.y Vertical location of the equation. Defaults to the top of the graph.
#' @param eq.size Size of the font of the equation.
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
#'   muck(600, 30, 2000, 0, 1, .01)
#' @export
muck <- function(y.obs, y.obs.er=0, y.target, y.target.er=0, F.contam=1, F.contam.er=0, MC=TRUE, its=1e4, filter=TRUE, roundby=1, decimals=3, visualise=TRUE, talk=TRUE, eq.x=0, eq.y=c(), eq.size=0.7, target.col="black", observed.col="blue", contamination.col="red", target.pch=20, observed.pch=18, contamination.pch=17, true.name="target", xlab="contamination (%)", ylab="F14C", ylim=c(), bty="l") {
  F.obs <- C14toF14C(y.obs, y.obs.er)
  F.target <- C14toF14C(y.target, y.target.er)

  # sanity checks
  if(all(F.contam > F.target[,1]) && all(F.obs[,1] < F.target[,1]))
    warning("The observed age (F.obs ", F.obs[,1], ") cannot be magicked into the younger target age (F.target ",
      F.target[,1], "), because the activity of the contamination (F.contam, ", F.contam, 
	  ") is higher than F.target.", call.=FALSE)
  if(all(F.contam < F.target[,1]) && all(F.obs[,1] > F.target[,1]))
    warning("The observed age (F.obs ", F.obs[,1], ") cannot be magicked into the older target age (F.target ",
      F.target[,1], "), because the activity of the contamination (F.contam, ", F.contam, 
	  ") is lower than F.target.", call.=FALSE)

  if(any(F.contam < 0))
    stop("F.contam cannot be smaller than 0", call.=FALSE)
  if(MC) {
    if(all(y.obs.er == 0)) # then waste no time sampling
	  F.observed.samples <- F.obs[1,1] else
        F.observed.samples <- rnorm(its, F.obs[,1], F.obs[,2])
    if(all(y.target.er == 0)) # then waste no time sampling
      F.target.samples <- F.target[1,1] else
        F.target.samples <- rnorm(its, F.target[,1], F.target[,2])
    if(all(F.contam.er == 0)) # then waste no time sampling
      F.contam.samples <- F.contam[1] else
        F.contam.samples <- rnorm(its, F.contam, F.contam.er)
    frac <- ((F.observed.samples - F.target.samples) / (F.contam.samples - F.target.samples))
    if(filter) {
      inside <- frac >= 0 & frac <= 1 # nothing outside [0,1]
	  frac <- frac[inside]
	  F.observed.samples <- F.observed.samples[inside]
	  F.target.samples <- F.target.samples[inside] 
	  F.contam.samples <- F.contam.samples[inside]
    }

    if(length(frac) < 10)
      stop("too few valid samples after filtering. Increase 'its'", call.=FALSE)
    perc.sd <- round(sd(100*frac), roundby)
    perc <- round(mean(100*frac), roundby)
    if(is.na(perc.sd)) # when all errors are 0, sd(0) becomes NA
      perc.sd <- 0
	
	F.contam <- mean(F.contam.samples)
	F.contam.er <- sd(F.contam.samples)
    if(is.na(F.contam.er)) # when all errors are 0, sd(0) becomes NA
      F.contam.er <- 0
	
  } else {
      # note: these calculations do NOT take into account uncertainties
      frac <- ((F.obs[,1] - F.target[,1]) / (F.contam - F.target[,1]))
      perc.sd <- 0
      perc <- 100*frac
    }

  if(visualise)
    if(length(y.obs) == 1) {
      plot_contamination(true.F=F.target[,1], true.er=F.target[,2], obs.F=F.obs[,1], obs.er=F.obs[,2], perc=perc, perc.er=perc.sd, contam.F=F.contam, contam.er=F.contam.er, ylim=ylim, xlab=xlab, true.col=target.col, observed.col=observed.col, contamination.col=contamination.col, true.pch=target.pch, true.name=true.name, observed.pch=observed.pch, contamination.pch=contamination.pch, ylab=ylab, bty=bty)
      txt <- c("contamination", " = (", round(F.obs[,1], decimals), "-", round(F.target[,1], decimals), ") / (", 
        round(F.contam, decimals), "-",  round(F.target[,1], decimals), ") = ", round(mean(frac), decimals))
      colours <- c(contamination.col, grey(.5), observed.col, grey(.5), target.col, grey(.5), contamination.col, grey(.5), target.col, grey(.5), contamination.col)
      xpos <- eq.x+(eq.size*cumsum(c(0, strwidth(txt[-length(txt)]))))
      if(length(eq.y) == 0)
        eq.y <- 1.01*max(c(F.obs, F.target[,1], F.contam))
      for(i in seq_along(txt))
        text(x = xpos[i], y = eq.y, labels = txt[i], col = colours[i], adj = c(0, 0), cex=eq.size)
	}  
 
  F.obs <- round(F.obs, decimals); F.target <- round(F.target, decimals); perc <- round(perc, decimals)

  if(talk)
    if(length(y.obs) == 1) {
      message("Observed age: ", y.obs, "+-", y.obs.er, " C14 BP (", F.obs[,1], "+-", F.obs[,2], " F14C)")
      message("Target age: ", y.target, " C14 BP (", F.target[,1], " F14C)")
      message("Calculation: (", F.obs[,1], "-", F.target[,1], ")/(", F.contam, " - ", F.target[,1], ") = ", round(mean(frac), decimals+2))
      if(MC)
        message("Contamination required: ", perc, "+-", perc.sd, "%" ) else
          message("Contamination required: ", perc, "%" )
    }
  if(MC)
    invisible(list(perc=perc, samples=cbind(F.observed.samples, F.target.samples, F.contam.samples, frac))) else
      invisible(cbind(perc))
}



#' @name push.normal
#' @title Add a normal distribution to a calibrated date
#' @description Push a date to younger or older ages by adding (or subtracting) a normal distribution (e.g. if a bone is assumed to have a lag or in-built age)
#' @details n random values will be sampled from the calibrated distribution, and a similar amount will be sampled from the normal distribution. The sampled values will then be added to or subtracted from each other to push the date to younger or older ages.
#' @return The resulting calibrated distribution and its hpd ranges, together with a plot of the pushed date with the normal distribution (and whether it is added or subtracted) as inset
#' @param y The radiocarbon age.
#' @param er The error of the radiocarbon age.
#' @param mean The mean of the normal or gamma distribution.
#' @param sdev The standard deviation of the normal distribution.
#' @param add The distribution can be added or subtracted. Adding results in ages being pushed to younger age distributions, and subtracting to older ones.
#' @param n The amount of random values to sample (from both the calibrated distribution and the gamma/normal distribution) to calculate the push. Defaults to \code{n=1e6}.
#' @param prob The probability for the hpd ranges. Defaults to \code{prob=0.95}.
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages. Defaults to \code{postbomb=FALSE}.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param cal.lim Calendar axis limits. Calculated automatically by default.
#' @param calib.col Colour of the calibrated distribution (defaults to semi-transparent light grey).
#' @param pushed.col Colour of the pushed distribution (defaults to semi-transparent blue).
#' @param inset Whether or not to plot an inset graph showing the shape of the normal/gamma distribution.
#' @param inset.col Colour of the normal/gamma distribution.
#' @param inset.loc Location of the inset graph.
#' @param inset.mar Margins of the inset graph.
#' @param inset.mgp Margin lines for the inset graph.
#' @examples
#'   push.normal(250, 25, 50, 10)
#' @export
push.normal <- function(y, er, mean, sdev, add=TRUE, n=1e6, prob=0.95, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, thiscurve=NULL, cc.dir=NULL, normal=TRUE, t.a=3, t.b=4, BCAD=FALSE, cal.lim=c(), calib.col=rgb(0,0,0,.25), pushed.col=rgb(0,0,1,.4), inset=TRUE, inset.col="darkgreen", inset.loc=c(0.6, 0.97, 0.6, 0.97), inset.mar=c(3, 0.5, 0.5, 0.5), inset.mgp=c(2,1,0)) {
  if(length(y) != 1 || length(er) != 1)
    stop("Please provide one value for both y and er")
  if(length(mean) != 1 || length(sdev) != 1)
    stop("Please provide one value for both mean and par2 (sdev or shape)")

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)
  
  shift <- rnorm(n, mean, sdev) 
  calib <- caldist(y, er, cc=cc, postbomb=postbomb, thiscurve=thiscurve, normalise=TRUE, BCAD=BCAD, cc.dir=cc.dir)
  rcalib <- r.calib(n, y, er, cc=cc, postbomb=postbomb, thiscurve=thiscurve, normal=normal, t.a=t.a, t.b=t.b, normalise=TRUE, BCAD=BCAD, rule=2, cc.dir=cc.dir)

  if(add) { # the date becomes younger
    shifted <- if(BCAD) rcalib + shift else rcalib - shift
   } else { # the date becomes older
       shifted <- if(BCAD) rcalib - shift else rcalib + shift
   }
  shifted <- density(shifted)
   
  if(length(cal.lim) == 0) {
    cal.lim <- range(calib[,1], shifted$x)
    if(!BCAD) cal.lim <- rev(cal.lim)
  }
  
  plot(0, type="n", xlim=cal.lim, ylim=c(0, 1.5), bty="l", ylab="", yaxt="s", xlab=ifelse(BCAD, "BCAD", "cal BP"))
  draw.dist(cbind(calib), dist.col=calib.col, dist.border=calib.col, y.pos=0, ex=-1)
  hpds <- draw.dist(cbind(shifted$x, shifted$y), prob=prob, y.pos=0, ex=-1)

  # inset graph
  if(inset) {
    op <- par(fig=inset.loc, new=TRUE, mar=inset.mar, mgp=inset.mgp, bty="n")
    xseq <- seq(mean-(3*sdev), mean+(3*sdev), length=200)
    if(add) xlim <- range(xseq) else xlim <- rev(range(xseq))
    plot(xseq, dnorm(xseq, mean, sdev), type="l", xlim=xlim, col=inset.col, xlab="", ylab="", yaxt="n", yaxs="r")
    end <- mean+sdev
    arrows(mean, 0, end, 0, col=inset.col, lwd=2, length=.05)
    op <- par(fig=c(0,1,0,1), mar=c(5,4,4,2), mgp=c(3,1,0))
  }
  print(hpds)
  invisible(list(shifted=shifted, hpds=hpds))
}



#' @name push.gamma
#' @title Add a gamma distribution to a calibrated date
#' @description Push a date to younger or older ages by adding (or subtracting) a gamma distribution (e.g. if a bone is assumed to have a lag or in-built age)
#' @details n random values will be sampled from the calibrated distribution, and a similar amount will be sampled from the gamma distribution. The sampled values will then be added to or subtracted from each other to push the date to younger or older ages.
#' @return The resulting calibrated distribution and its hpd ranges, together with a plot of the pushed date with the gamma distribution (and whether it is added or subtracted) as inset
#' @param y The radiocarbon age
#' @param er The error of the radiocarbon age
#' @param mean The mean of the gamma distribution
#' @param shape The shape of the gamma distribution. If setting this to shape=1, it becomes an exponential distribution.
#' @param add The distribution can be added or subtracted. Adding results in ages being pushed to younger age distributions, and subtracting to older ones.
#' @param n The amount of random values to sample (from both the calibrated distribution and the gamma distribution) to calculate the push. Defaults to \code{n=1e6}.
#' @param prob The probability for the hpd ranges. Defaults to \code{prob=0.95}.
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages. Defaults to \code{postbomb=FALSE}.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param cal.lim Calendar axis limits. Calculated automatically by default.
#' @param calib.col Colour of the calibrated distribution (defaults to semi-transparent light grey).
#' @param pushed.col Colour of the pushed distribution (defaults to semi-transparent blue).
#' @param inset Whether or not to plot an inset graph showing the shape of the normal/gamma distribution.
#' @param inset.col Colour of the normal/gamma distribution.
#' @param inset.loc Location of the inset graph.
#' @param inset.mar Margins of the inset graph.
#' @param inset.mgp Margin lines for the inset graph.
#' @examples
#'   push.gamma(250, 25, 50, 2, add=FALSE) # subtract a gamma distribution
#' @export
push.gamma <- function(y, er, mean, shape, add=TRUE, n=1e6, prob=0.95, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, thiscurve=NULL, cc.dir=NULL, normal=TRUE, t.a=3, t.b=4, BCAD=FALSE, cal.lim=c(), calib.col=rgb(0,0,0,.25), pushed.col=rgb(0,0,1,.4), inset=TRUE, inset.col="darkgreen", inset.loc=c(0.6, 0.97, 0.6, 0.97), inset.mar=c(3, 0.5, 0.5, 0.5), inset.mgp=c(2,1,0)) {
  if(length(y) != 1 || length(er) != 1)
    stop("Please provide one value for both y and er")
  if(length(mean) != 1 || length(shape) != 1)
    stop("Please provide one value for both mean and par2 (sdev or shape)")

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)
  
  shift <- rgamma(n, shape, shape/mean)
  calib <- caldist(y, er, cc=cc, postbomb=postbomb, thiscurve=thiscurve, normalise=TRUE, BCAD=BCAD, cc.dir=cc.dir)
  rcalib <- r.calib(n, y, er, cc=cc, postbomb=postbomb, thiscurve=thiscurve, normal=normal, t.a=t.a, t.b=t.b, normalise=TRUE, BCAD=BCAD, rule=2, cc.dir=cc.dir)

  if(add) { # the date becomes younger
    shifted <- if(BCAD) rcalib + shift else rcalib - shift
   } else { # the date becomes older
       shifted <- if(BCAD) rcalib - shift else rcalib + shift
   }
  shifted <- density(shifted)
   
  calpol <- cbind(c(calib[,1], rev(calib[,1])), c(calib[,2]/max(calib[,2]), rep(0, nrow(calib))))
  shiftpol <- cbind(c(min(shifted$x), shifted$x, max(shifted$x)), c(0, shifted$y/max(shifted$y), 0))

  if(length(cal.lim) == 0) {
    cal.lim <- range(calib[,1], shifted$x)
    if(!BCAD) cal.lim <- rev(cal.lim)
  }
  
  plot(0, type="n", xlim=cal.lim, ylim=c(0, 1.5), bty="l", ylab="", yaxt="s", xlab=ifelse(BCAD, "BCAD", "cal BP"))
  draw.dist(cbind(calib), dist.col=calib.col, dist.border=calib.col, y.pos=0, ex=-1)
  hpds <- draw.dist(cbind(shifted$x, shifted$y), prob=prob, y.pos=0, ex=-1)
 
  # inset graph
  if(inset) {
    op <- par(fig=inset.loc, new=TRUE, mar=inset.mar, mgp=inset.mgp, bty="n")
    xseq <- seq(0, mean*(1+4/sqrt(shape)), length=200)
    if(add) xlim <- range(xseq) else xlim <- rev(range(xseq))
    plot(xseq, dgamma(xseq, shape, shape/mean), type="l", xlim=xlim, col=inset.col, xlab="", ylab="", yaxt="n", yaxs="r")
    end <- mean+(mean/3)
    arrows(mean, 0, end, 0, col=inset.col, lwd=2, length=.05)
    op <- par(fig=c(0,1,0,1), mar=c(5,4,4,2), mgp=c(3,1,0))
  }
  print(hpds)
  invisible(list(shifted=shifted, hpds=hpds))
}


