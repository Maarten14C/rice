#' @name caldist
#' @title Calculate calibrated distribution
#' @description Calculate the calibrated distribution of a radiocarbon date.
#' @return The probability distribution(s) as two columns: cal BP ages and their associated probabilities
#' @param y Uncalibrated radiocarbon age
#' @param er Lab error of the radiocarbon age
#' @param cc Calibration curve to use. Defaults to IntCal20 (cc=1), can be Marine20 (cc=2), SHCal20 (cc=3), or if postbomb=TRUE, NH1 (cc=1), NH2 (cc=2), NH3 (cc=3), SH1-2 (cc=4) or SH3 (cc=5).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param bombalert Stop if a date is overly close to the younger limit of the IntCal curve. Defaults to \code{bombalert=TRUE}. This error can be avoided by either providing a postbomb curve (e.g., \code{postbomb=1}) or typing \code{bombalert=FALSE} (in this case, part of the date will be truncated).
#' @param glue Glue postbomb and prebomb curves together. Defaults to 0 (none), can be 1 (IntCal20 + NH1), 2 (IntCal20 + NH2), 3 (IntCal20 + NH3), 4 (SHCal20 + SH1-2) or 5 (SHCal20 + SH3). Note that this will override the value of cc.
#' @param deltaR Age offset (e.g. for marine samples). This assumes that the radiocarbon age is provided as 14C BP (not F14C or pMC).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param is.F Set this to TRUE if the provided age and error are in the F14C timescale.
#' @param is.pMC Set this to TRUE if the provided age and error are in the pMC timescale.
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=TRUE}, so not using the C14 timescale (will be more accurate especially for dates with larger errors, e.g., older ones).
#' @param thiscurve As an alternative to providing cc and/or postbomb/glue, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param yrsteps Steps to use for interpolation. Defaults to the cal BP steps in the calibration curve
#' @param cc.resample The IntCal20 curves have different densities (every year between 0 and 5 kcal BP, then every 5 yr up to 15 kcal BP, then every 10 yr up to 25 kcal BP, and then every 20 yr up to 55 kcal BP). If calibrated ages span these density ranges, their drawn heights can differ, as can their total areas (which should ideally all sum to the same size). To account for this, resample to a constant time-span, using, e.g., \code{cc.resample=5} for 5-yr timespans.
#' @param dist.res As an alternative to yrsteps, provide the amount of 'bins' in the distribution.
#' @param pb.steps Yearly steps for postbomb curves. Defaults to 20 steps per year, \code{pb.steps=0.05}.
#' @param cc0.res Length of 'curve' when cc=0 (no calibration curve). Defaults to 5000, in order to provide enough points for detailed distributions.  
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=1e-6}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param normalise Sum the entire calibrated distribution to 1. Defaults to \code{normalise=TRUE}.
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param rule Which extrapolation rule to use. Defaults to \code{rule=1} which returns NAs.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param col.names Names for the output columns. Defaults to calBP/BCAD and probs, respectively (depending on the value of BCAD).
#' @examples
#' calib <- caldist(130,10, bombalert=FALSE)
#' plot(calib, type="l")
#' postbomb <- caldist(-3030, 20, postbomb=1, BCAD=TRUE)
#' @export
caldist <- function(y, er, cc=1, postbomb=FALSE, bombalert=TRUE, glue=0, deltaR=0, deltaSTD=0, is.F=FALSE, is.pMC=FALSE, as.F=TRUE, thiscurve=NULL, yrsteps=FALSE, cc.resample=FALSE, dist.res=200, pb.steps=0.05, cc0.res=5e3, threshold=1e-3, normal=TRUE, t.a=3, t.b=4, normalise=TRUE, BCAD=FALSE, rule=1, cc.dir=NULL, col.names=NULL) {

  if(is.F && is.pMC)
    stop("cannot have both is.F=TRUE and is.pMC=TRUE")

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  if(length(thiscurve) > 0)
    this.cc <- thiscurve else {
      if(cc == 0) { # no ccurve needed
        xseq <- seq(y-4*er, y+4*er, length=cc0.res)
        this.cc <- cbind(xseq, xseq, rep(0, length(xseq)))
      } else {
          if(glue>0) {
            if(glue %in% 1:3)
              this.cc <- rintcal::glue.ccurves(1, postbomb=glue, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                if(glue %in% 4:5)			
                  this.cc <- rintcal::glue.ccurves(3, postbomb=glue, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                    stop("please provide an integer for glue between 0 and 5")
		  } else {
		  
            this.cc <- rintcal::ccurve(cc, postbomb=postbomb, cc.dir=cc.dir, 
              resample=cc.resample, as.F=is.F, as.pMC=is.pMC)
         
            # check if any dates lie at the younger edge of this.cc
            young <- FALSE
            if(is.F || is.pMC) {
              if(max(y+(3*er)) > max(this.cc[,2]))
                young <- TRUE
            } else 
                if(min(y-(3*er)) < min(this.cc[,2]))
                  young <- TRUE 

            if(young) {
              if(postbomb == FALSE) {
                if(bombalert)
                  stop("please provide a postbomb curve") else { 
                    # no postbomb has been defined, but we need to add one 
                    if(cc==1) # then assume we need NH1
                      this.cc <- rintcal::glue.ccurves(cc, postbomb=1, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                      if(cc==3) # then assume we need SH1-2
                        this.cc <- rintcal::glue.ccurves(cc, postbomb=4, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                          stop("please provide a postbomb curve, e.g. postbomb=1")
                } 
              } else
                  this.cc <- rintcal::glue.ccurves(cc, postbomb=postbomb, cc.dir, as.F=is.F, as.pMC=is.pMC)
            }
          }
		}  
    }

  if(postbomb) {
    xseq <- seq(min(this.cc[,1]), max(this.cc[,1]), by=pb.steps)
    ccmu <- approx(this.cc[,1], this.cc[,2], xseq)$y
    ccsd <- approx(this.cc[,1], this.cc[,3], xseq)$y
    this.cc <- cbind(xseq, ccmu, ccsd)
  }

  if(as.F && !is.F) { # y, er and cc are in C14 timescale, but need to be in F
    this.cc <- cbind(this.cc[,1], C14toF14C(this.cc[,2], this.cc[,3]))
    asF <- as.numeric(C14toF14C(y, er))
    y <- asF[1]; er <- asF[2]
  }

  # calibrate; find how far age (measurement) is from cc[,2] of calibration curve
  if(normal)
    cal <- cbind(this.cc[,1], dnorm(this.cc[,2], y, sqrt(er^2+this.cc[,3]^2))) else
      cal <- cbind(this.cc[,1], ((t.b + ((y-this.cc[,2])^2) / (2*(this.cc[,3]^2 + er^2))) ^ (-1*(t.a+0.5)))/sqrt(this.cc[,3]^2 + er^2) ) # corrected 13 Jan 2026

  # interpolate and normalise calibrated distribution to 1
  if(postbomb)
    if(!yrsteps)
      yrsteps <- 0.05 # enough detail to enable calculation of hpd ranges also for postbomb dates
  if(yrsteps)
    yrseq <- seq(min(cal[,1], na.rm=TRUE), max(cal[,1], na.rm=TRUE), by=yrsteps) else
      yrseq <- cal[,1]
 #     yrsteps <- seq(min(cal[,1]), max(cal[,1]), length=dist.res)
  cal <- approx(cal[,1], cal[,2], yrseq, rule=rule)
  # cal <- cbind(cal$x, cal$y/sum(cal$y)) # normalise
  cal <- cbind(cal$x, cal$y)

  if(normalise)
    cal[,2] <- cal[,2]/sum(cal[,2])
  # remove years with very small probabilities on the extremes of the distribution
  above <- which(cal[,2] >= (threshold * max(cal[,2]))) # relative to its peak
  if(length(above)>2)
    cal <- cal[min(above):max(above),] # now does not necessarily sum to exactly 1 any more
  if(normalise) # ... so normalise again if asked for
    cal[,2] <- cal[,2]/sum(cal[,2])

  if(is.null(col.names)) {
    colnames(cal) <- c("cal BP", "prob")
    if(BCAD) {
      cal[,1] <- calBPtoBCAD(cal[,1])
      colnames(cal)[1] <- "cal BC/AD"
    }
  } else
      if(!is.na(col.names[1]))
        if(length(col.names) == 2)
          colnames(cal) <- col.names

  return(cal)
}



#' @name point.estimates
#' @title Calculate a point estimate
#' @description Calculate a point estimate of a calibrated distribution - either the weighted mean, the median or the mode (maximum). Note that point estimates often tend to be very poor representations of entire calibrated distributions, so please be careful and do not reduce entire calibrated distributions to just 1 point value.
#' @return The chosen point estimates
#' @param calib The calibrated distribution, as returned from caldist()
#' @param wmean Report the weighted mean (defaults to TRUE)
#' @param median Report the median (defaults to TRUE)
#' @param mode Report the mode, which is the year with the maximum probability (defaults to TRUE)
#' @param midpoint Report the midpoint of the hpd range(s)
#' @param prob probability range for the hpd range(s)
#' @param rounded Rounding for reported probabilities. Defaults to 1 decimal.
#' @param every Yearly precision (defaults to \code{every=1}).
#' @examples
#' point.estimates(caldist(130,20, bombalert=FALSE))
#' plot(tmp <- caldist(2450,50), type='l')
#' abline(v=point.estimates(tmp), col=1:4)
#' @export
point.estimates <- function(calib, wmean=TRUE, median=TRUE, mode=TRUE, midpoint=TRUE, prob=.95, rounded=1, every=1) {
  to.report <- c()
  name <- c()
  calib[,2] <- calib[,2] / sum(calib[,2])

  if(wmean) {
    wmean <- weighted.mean(calib[,1], calib[,2])
    to.report <- c(to.report, wmean)
    name <- c(name, "weighted mean")
  }
  if(median) {
    median <- approx(cumsum(calib[,2]), calib[,1], 0.5, rule=2)$y
    to.report <- c(to.report, median)
    name <- c(name, "median")
  }
  if(mode) {
    mode <- calib[which(calib[,2] == max(calib[,2]))[1],1]
    to.report <- c(to.report, mode)
    name <- c(name, "mode")
  }
  if(midpoint) {
    midpoint <- range(hpd(calib, prob, prob.round=rounded, every=every)[,1:2])
    midpoint <- midpoint[1] + (midpoint[2]-midpoint[1])/2
    to.report <- c(to.report, midpoint)
    name <- c(name, "midpoint")
  }

  names(to.report) <- name
  return(round(to.report, rounded))
}



#' @name hpd
#' @title Calculate highest posterior density
#' @description Calculate highest posterior density ranges of a calibrated distribution
#' @return The highest posterior density ranges, as three columns: from age, to age, and the corresponding percentage(s) of the range(s)
#' @param calib The calibrated distribution, as returned from caldist()
#' @param prob Probability range which should be calculated. Default \code{prob=0.95}.
#' @param return.raw The raw data to calculate hpds can be returned, e.g. to draw polygons of the calibrated distributions. Defaults to \code{return.raw=FALSE}.
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param ka Whether to report results in years (default) or as ka
#' @param age.round Rounding for ages. Defaults to 0 decimals.
#' @param prob.round Rounding for reported probabilities. Defaults to 1 decimal.
#' @param every Yearly precision (defaults to 0.1, as a compromise between speed and accuracy).
#' @param bins The number of bins required. Any distribution with fewer bins gets recalculated using 100 narrower bins.
#' @examples
#' hpd(caldist(130,20, bombalert=FALSE))
#' plot(tmp <- caldist(2450,50), type='l')
#' myhpds <- hpd(tmp)
#' abline(v=unlist(myhpds[,1:2]), col=4)
#' @export
hpd <- function(calib, prob=0.95, return.raw=FALSE, BCAD=FALSE, ka=FALSE, age.round=0, prob.round=1, every=0.1, bins=20) {

  # re-interpolate to desired precision
  if(ka) {
    every <- every/1e3
    age.round <- age.round+3
  }

  calib <- calib[order(calib[,1]),] # ensure calendar ages are in increasing order
  steppies <- diff(unique(calib[,1]))
  tmpcalib <- approx(calib[,1], calib[,2], seq(min(calib[,1]), max(calib[,1]), by=every))
  if(length(tmpcalib$x) >= bins) # then the distribution is not very narrow
    calib <- tmpcalib else
      calib <- approx(calib[,1], calib[,2], seq(min(calib[,1]), max(calib[,1]), length=100), rule=2)
  calib <- cbind(calib$x, calib$y/sum(calib$y)) # normalise probs to 1

  # rank the calibrated ages according to their probabilities
  o <- order(calib[,2], decreasing=TRUE) # decreasing probs
  calib <- cbind(calib[o,], cumsum(calib[o,2]) / sum(calib[,2]))
  calib <- cbind(calib, calib[,3] <= prob) # yr, probs, cumprobs, within/outside ranges (T/F)
  calib <- calib[order(calib[,1], decreasing=TRUE),] # put ages in order again

  # find the outer ages of the calibrated ranges
  if(BCAD) {
    from <- sort(calib[which(diff(c(calib[,4], 0)) == -1), 1])
    to <- sort(calib[which(diff(c(0, calib[,4])) == 1), 1])

    # swap order when the range is in AD or crosses BC/AD
    swap <- which( (pmax(from, to) > 0) & (from > to) )
    if(length(swap) > 0) {
      tmp <- from[swap]
      from[swap] <- to[swap]
      to[swap] <- tmp
    }
  } else {
    from <- sort(calib[which(diff(c(0, calib[,4])) == 1), 1])
    to <- sort(calib[which(diff(c(calib[,4], 0)) == -1), 1])
  }

  # find the year and probability within each range (as %)
  perc <- 0
  for(i in 1:length(from))
    if(from[i] == to[i])
      perc[i] <- calib[which(calib[,1] == from[i]),2] else {
        fromto <- calib[which(calib[,1] == from[i])[1] : which(calib[,1]== to[i])[1],1:2]
        perc[i] <- round(100*sum(fromto[,2]), prob.round)
      }

  mindiff <- min(abs(diff(calib[,1])))
  if(mindiff < 0.01) # very small values require different rounding
    age.round <- 1+nchar(gsub("[1-9].*", "", sub("0\\.", "",
      format(mindiff, scientific = FALSE))))

  hpds <- data.frame(from=round(from, age.round), to=round(to, age.round), 
    perc=round(perc, prob.round))

  if(return.raw)
    return(list(calib=calib[,-3], hpds=hpds)) else
      return(data.frame(hpds))
}




#' @name age.range
#' @title Calculate age ranges
#' @description Calculate the quantile age ranges of a calibrated distribution
#' @return The highest posterior density ranges, as three columns: from age, to age, and the corresponding percentage(s) of the range(s)
#' @param calib The calibrated distribution, as returned from caldist()
#' @param prob Probability range which should be calculated. Default \code{prob=0.95}.
#' @param roundby Rounding. Defaults to 0 decimals.
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @examples
#' age.range(caldist(130,20, bombalert=FALSE))
#' @export
age.range <- function(calib, prob=0.95, roundby=0, BCAD=FALSE) {
  if(NCOL(calib) == 1) { # then it's NOT a 'calib' and we assume it's a vector
    rng <- quantile(calib, c((1-prob)/2, 1 - (1-prob)/2))
  }	else {
    cdf <- cumsum(calib[,2] / sum(calib[,2]))
    rng <- approx(cdf, calib[,1],
      c((1-prob)/2, 1 - (1-prob)/2))$y
    }
  rng <- rev(rng)
  if(BCAD && max(rng) > 0)
    rev(rng)
  return(round(rng, roundby))
}


#  find the calibrated probability of a calendar age for a 14C date
#' @name l.calib
#' @title Find the calibrated probability of a calendar age for a 14C date.
#' @description Find the calibrated probability of a cal BP age for a radiocarbon date. Can handle either multiple calendar ages for a single radiocarbon date, or a single calendar age for multiple radiocarbon dates.
#' @details The function cannot deal with multiple calibration curves if multiple calendar years or radiocarbon dates are entered.
#' @return The calibrated probability of a calendar age for a 14C age
#' @param x The cal BP year.
#' @param y The radiocarbon date's mean.
#' @param er The radiocarbon date's lab error.
#' @param cc calibration curve for the radiocarbon date(s) (see the \code{rintcal} package).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param deltaR Age offset (e.g. for marine samples). This assumes that the radiocarbon age is provided as 14C BP (not F14C or pMC).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param is.F Use this if the provided date is in the F14C timescale.
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @author Maarten Blaauw
#' @examples
#'   l.calib(100, 130, 20)
#'   l.calib(100:110, 130, 20) # multiple calendar ages of a single date
#'   l.calib(100, c(130,150), c(15,20)) # multiple radiocarbon ages and a single calendar age
#'   plot(0:300, l.calib(0:300, 130, 20), type='l')
#' @export
l.calib <- function(x, y, er, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, thiscurve=c(), cc.dir=c(), normal=TRUE, as.F=FALSE, is.F=FALSE, t.a=3, t.b=4) {
  
  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  if(length(x) == 1) {
    if(length(y) != length(er))
      stop("check that y has as many entries as er") 		  
  } else
    if(length(y) > 1 || length(er)>1)
      stop("cannot deal with multiple entries for both x and y+er") 
  
  if(cc == 0) { # June 2025
    if(normal)
      prob <- dnorm(x, y, er) else 
        prob <- ( (t.b + ((x - y)^2) / (2 * er^2)) ^ (-1 * (t.a + 0.5))) /er
    prob[is.na(prob)] <- 0
    return(prob)	
  } else {
    if(is.F) {
      mu <- calBPtoF14C(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)
    } else
        if(as.F) {
          mu <- calBPtoF14C(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)
          tmp <- C14toF14C(y, er)
          y <- tmp[,1]; er <- tmp[,2]
        } else
          mu <- calBPtoC14(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)
    if(normal)
      prob <- dnorm(y, mu[,1], sqrt(mu[,2]^2 + er^2)) else
        prob <- ((t.b + ((y-mu[,1])^2) / (2*(sqrt(er^2+mu[,2]^2)^2))) ^ (-1*(t.a+0.5))) / sqrt(er^2+mu[,2]^2)
    prob[is.na(prob)] <- 0
    return(prob)
  }
}



#' @name r.calib
#' @title return a random calendar age from a calibrated distribution
#' @description Calculate the cumulative calibrated distribution, then sample n random uniform values between 0 and 1 and find the corresponding calendar ages through interpolation. Calendar ages with higher calibrated probabilities will be proportionally more likely to be sampled.
#' @return n randomly sampled calendar ages
#' @param n The number of calendar ages to sample
#' @param y Uncalibrated radiocarbon age
#' @param er Lab error of the radiocarbon age
#' @param cc Calibration curve to use. Defaults to IntCal20 (\code{cc=1}).
#' @param postbomb Whether or not to use a postbomb curve. Required for negative radiocarbon ages.
#' @param bombalert Warn if a date is close to the lower limit of the calibration curve. Defaults to \code{postbomb=TRUE}.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param is.F Use this if the provided date is in the F14C timescale.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param yrsteps Steps to use for interpolation. Defaults to the cal BP steps in the calibration curve
#' @param cc.resample The IntCal20 curves have different densities (every year between 0 and 5 kcal BP, then every 5 yr up to 15 kcal BP, then every 10 yr up to 25 kcal BP, and then every 20 yr up to 55 kcal BP). If calibrated ages span these density ranges, their drawn heights can differ, as can their total areas (which should ideally all sum to the same size). To account for this, resample to a constant time-span, using, e.g., \code{cc.resample=5} for 5-yr timespans.
#' @param dist.res As an alternative to yrsteps, provide the amount of 'bins' in the distribution
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=0}.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param normalise Sum the entire calibrated distribution to 1. Defaults to \code{normalise=TRUE}.
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param rule Which extrapolation rule to use. Defaults to \code{rule=1} which returns NAs.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param seed For reproducibility, a seed can be set (e.g., \code{seed=123}). Defaults to NA, no seed set.
#' @author Maarten Blaauw
#' @examples
#'   r.calib(10,130,20, bombalert=FALSE) # 10 random cal BP ages
#'   plot(density(r.calib(1e6, 2450, 20)))
#' @export
r.calib <- function(n, y, er, cc=1, postbomb=FALSE, bombalert=TRUE, deltaR=0, deltaSTD=0, as.F=FALSE, is.F=FALSE, thiscurve=NULL, yrsteps=FALSE, cc.resample=FALSE, dist.res=200, threshold=0, normal=TRUE, t.a=3, t.b=4, normalise=TRUE, BCAD=FALSE, rule=2, cc.dir=NULL, seed=NA) {
  if(length(n) == 0 || n<1)
    stop("n needs to be a value >0")
  if(!length(y) == 1 || !length(er) == 1)
    stop("I can only handle one date at a time")
  
  calib <- caldist(y, er, cc=cc, postbomb=postbomb, bombalert=bombalert, deltaR=deltaR, deltaSTD=deltaSTD, as.F=as.F, is.F=is.F, thiscurve=thiscurve, yrsteps=yrsteps, normalise=normalise, BCAD=BCAD, rule=rule, cc.dir=cc.dir)
  if(!is.na(seed))
    if(is.numeric(seed))  
      set.seed(seed) else
        message("seed has to be numeric")
  rx <- approx(cumsum(calib[,2])/sum(calib[,2]), calib[,1], runif(n), rule=2)$y
  return(rx)
}



#' @name younger
#' @title Find the probability of a calibrated date being of a certain age or younger than it
#' @description Find the probability that a sample is of a certain calendar age x or younger than it, by calculating the proportion of the calibrated distribution up to and including x (i.e., summing the calibrated distribution up to year x).
#' @details The function can only deal with one date at a time.
#' @return The probability of a date being of a certain calendar age or younger than it.
#' @param x The year of interest, in cal BP by default.
#' @param y The radiocarbon date's mean.
#' @param er The radiocarbon date's lab error.
#' @param cc calibration curve for the radiocarbon date(s) (see the \code{rintcal} package).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param bombalert Warn if a date is close to the lower limit of the calibration curve. Defaults to \code{postbomb=TRUE}.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param is.F Use this if the provided date is in the F14C timescale.
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=0}.
#' @author Maarten Blaauw
#' @examples
#' younger(2800, 2450, 20)
#' younger(2400, 2450, 20)
#' calibrate(160, 20, BCAD=TRUE)
#' younger(1750, 160, 20, BCAD=TRUE)
#' @export
younger <- function(x, y, er, cc=1, postbomb=FALSE, bombalert=TRUE, deltaR=0, deltaSTD=0, normal=TRUE, as.F=FALSE, is.F=FALSE, t.a=3, t.b=4, BCAD=FALSE, threshold=0) {
  if(length(y)>1 || length(er) >1)
    stop("I can only deal with one date at a time")
  
  cal <- caldist(y, er, cc, postbomb=postbomb, bombalert=bombalert, deltaR=deltaR, deltaSTD=deltaSTD, normal=normal, t.a=t.a, t.b=t.b, as.F=as.F, is.F=is.F, BCAD=BCAD, threshold=threshold)

  dx <- diff(cal[,1])
  cdf <- c(0, cumsum((cal[-nrow(cal),2] + cal[-1,2]) / 2 * dx)) # trapezium integration
  cdf <- cdf / max(cdf)

  return(approx(cal[,1], cdf, x, rule=2)$y) # cumulative prob
}



#' @name older
#' @title Find the probability of a calibrated date being older than a certain age
#' @description Find the probability of a calibrated date being older than an age x.
#' @description Find the probability that a sample is older than a certain calendar age x, by calculating the proportion of the calibrated distribution 'after' x (i.e., 1 - the summed calibrated distribution up to year x).
#' @details The function can only deal with one date at a time.
#' @return The probability of a date being older than a certain calendar age.
#' @param x The year of interest, in cal BP by default.
#' @param y The radiocarbon date's mean.
#' @param er The radiocarbon date's lab error.
#' @param cc calibration curve for the radiocarbon date(s) (see the \code{rintcal} package).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param bombalert Warn if a date is close to the lower limit of the calibration curve. Defaults to \code{postbomb=TRUE}.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param is.F Use this if the provided date is in the F14C timescale.
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=0}.
#' @author Maarten Blaauw
#' @examples
#' older(2800, 2450, 20)
#' older(2400, 2450, 20)
#' calibrate(160, 20, BCAD=TRUE)
#' older(1750, 160, 20, BCAD=TRUE)
#' @export
older <- function(x, y, er, cc=1, postbomb=FALSE, bombalert=TRUE, deltaR=0, deltaSTD=0, normal=TRUE, as.F=FALSE, is.F=FALSE, t.a=3, t.b=4, BCAD=FALSE, threshold=0) {
  return(1 - younger(x, y, er, cc, postbomb=postbomb, bombalert=bombalert, deltaR, deltaSTD, normal, as.F=as.F, is.F=is.F, t.a=t.a, t.b=t.b, BCAD=BCAD, threshold=threshold))
}



#' @name p.range
#' @title Probability of a date lying within a cal BP range
#' @description Find the probability of a calibrated date lying within an age range
#' @details The function can only deal with one date at a time.
#' @return The probability of a date lying within a certain calendar age range.
#' @param x1 The start the range of interest.
#' @param x2 The end of the range of interest.
#' @param y The radiocarbon date's mean.
#' @param er The radiocarbon date's lab error.
#' @param cc calibration curve for the radiocarbon date(s) (see the \code{rintcal} package).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param bombalert Warn if a date is close to the lower limit of the calibration curve. Defaults to \code{postbomb=TRUE}.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param is.F Use this if the provided date is in the F14C timescale.
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=0}.
#' @author Maarten Blaauw
#' @examples
#' p.range(2800, 2400, 2450, 20)
#' @export
p.range <- function(x1, x2, y, er, cc=1, postbomb=FALSE, bombalert=TRUE, deltaR=0, deltaSTD=0, normal=TRUE, as.F=FALSE, is.F=FALSE, t.a=3, t.b=4, BCAD=FALSE, threshold=0) {
  if(length(y)>1 || length(er) >1)
    stop("I can only deal with one date at a time")
  cal <- caldist(y, er, cc, postbomb=postbomb, bombalert=bombalert, deltaR=deltaR, deltaSTD=deltaSTD, normal=normal, t.a=t.a, t.b=t.b, as.F=as.F, is.F=is.F, BCAD=BCAD, threshold=threshold)

  dx <- diff(cal[,1])
  cdf <- c(0, cumsum((cal[-nrow(cal),2] + cal[-1,2]) / 2 * dx)) # trapezium integration
  cdf <- cdf / max(cdf)

  prob <- approx(cal[,1], cdf, sort(c(x1, x2)), rule=2)$y
  #prob <- approx(cal[,1], cumsum(cal[,2])/sum(cal[,2]), sort(c(x1, x2)), rule=2)$y
  return(max(prob)-min(prob))
}



#' @name calib.t
#' @title Comparison dates calibrated using both the t distribution (Christen and Perez 2009) and the normal distribution.
#' @description Visualise how a date calibrates using the t distribution and the normal distribution.
#' @details Radiocarbon and other dates are usually modelled using the normal distribution (red curve). The t approach (grey distribution) however allows for wider tails and thus tends to better accommodate outlying dates. This distribution requires two parameters, called 'a' and 'b'.
#' @param y The reported mean of the date.
#' @param er The reported error of the date.
#' @param t.a Value for the t parameter \code{a}.
#' @param t.b Value for the t parameter \code{b}.
#' @param cc calibration curve for the radiocarbon date(s) (see the \code{rintcal} package).
#' @param postbomb Which postbomb curve to use for negative 14C dates.
#' @param deltaR Age offset (e.g. for marine samples).
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation).
#' @param as.F Whether or not to calculate ages in the F14C timescale. Defaults to \code{as.F=FALSE}, which uses the C14 timescale.
#' @param is.F Use this if the provided date is in the F14C timescale.
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}.
#' @param cal.rev Reverse the calendar age axis. Defaults to TRUE
#' @param cc.dir Directory where the calibration curves for C14 dates \code{cc} are allocated. By default \code{cc.dir=c()}.
#' Use \code{cc.dir="."} to choose current working directory. Use \code{cc.dir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param normal.col Colour of the normal curve
#' @param normal.lwd Line width of the normal curve
#' @param t.col Colour of the t histogram
#' @param t.border Colour of the border of the t histogram
#' @param xlim x axis limits
#' @param ylim y axis limits
#' @author Maarten Blaauw
#' @examples
#' calib.t()
#'
#' @export
calib.t <- function(y=2450, er=50, t.a=3, t.b=4, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, as.F=FALSE, is.F=FALSE, BCAD=FALSE, cal.rev=TRUE, cc.dir=c(), normal.col="red", normal.lwd=1.5, t.col=rgb(0,0,0,0.25), t.border=rgb(0,0,0,0,0.25), xlim=c(), ylim=c()) {

  y <- y - deltaR
  er <- sqrt(er^2 + deltaSTD^2)

  normalcal <- caldist(y, er, cc, postbomb, BCAD=BCAD, cc.dir=cc.dir, as.F=as.F, is.F=is.F, normal=TRUE)
  tcal <- caldist(y, er, cc, postbomb, BCAD=BCAD, as.F=as.F, is.F=is.F, t.a=t.a, t.b=t.b, cc.dir=cc.dir, normal=FALSE)
  tpol <- cbind(c(tcal[1,1], tcal[,1], tcal[nrow(tcal),1]), c(0, tcal[,2], 0))

  xlab <- ifelse(BCAD, "BC/AD", "cal BP")
  if(length(xlim) == 0)
    xlim <- range(c(tpol[,1], normalcal[, 1]))
  if(cal.rev)
    xlim <- rev(xlim)
  if(BCAD)
    xlim <- rev(xlim)
  if(length(ylim) == 0)
    ylim <- c(0, max(tcal[,2], normalcal[, 2]))
  plot(normalcal, type="l", xlab=xlab, xlim=xlim, ylab="", ylim=ylim, col=normal.col, lwd=normal.lwd)
  polygon(tpol, col=t.col, border=t.border)
  legend("topleft", "Gaussian", text.col=normal.col, bty="n")
  legend("topright", paste("t (a=", t.a, ", b=", t.b, ")", sep=""), bty="n", text.col=t.col)
}



#' @name smooth.curve
#' @title Smooth a calibration curve
#' @description Smooth a calibration curve over a time window of a specified width. This to accommodate material that has accumulated over a certain assumed time, e.g. a cm of peat over say 30 years.
#' @details The smoothing is done by calculating the mean C14 age and error of a moving window (moving along with the cal BP steps of the calibration curve). Something similar is done in the online calibration software CALIB.
#' @param smooth The window width of the smoothing. Defaults to \code{smooth=30}.
#' @param cc The calibration curve to smooth. Calibration curve for 14C dates: 'cc=1' for IntCal20 (northern hemisphere terrestrial), 'cc=2' for Marine20 (marine), 'cc=3' for SHCal20 (southern hemisphere terrestrial). Alternatively, one can also write, e.g., "IntCal20", "Marine13". One can also make a custom-built calibration curve, e.g. using 'mix.ccurves()', and load this using 'cc=4'. In this case, it is recommended to place the custom calibration curve in its own directory, using 'cc.dir' (see below).
#' @param postbomb Use 'postbomb=TRUE' to get a postbomb calibration curve (default 'postbomb=FALSE'). For monthly data, type e.g. 'ccurve("sh1-2_monthly")'
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., 'cc.dir="ccurves"'.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). Defaults to c().
#' @param resample The IntCal curves come at a range of 'bin sizes'; every year from 0 to 5 kcal BP, then every 5 yr until 15 kcal BP, then every 10 yr until 25 kcal BP, and every 20 year thereafter. The curves can be resampled to constant bin sizes, e.g. 'resample=5'. Defaults to FALSE.
#' @param name The filename of the curve, if it is being saved. Defaults to \code{name="smoothed.csv"}.
#' @param save Whether or not to save the curve to cc.dir. Defaults to \code{save=FALSE}.
#' @param sep Separator between fields if the file is saved (tab by default, \code{sep="\t"}).
#' @author Maarten Blaauw
#' @examples
#'  mycurve <- smooth.ccurve(smooth=50)
#'  calibrate(2450,20, thiscurve=mycurve)
#' @export
smooth.ccurve <- function(smooth=30, cc=1, postbomb=FALSE, cc.dir=c(), thiscurve=c(), resample=0, name="smoothed.csv", save=FALSE, sep="\t") {
  Cc <- rintcal::ccurve(cc, postbomb=postbomb, cc.dir=cc.dir, resample=resample)
  if(length(thiscurve) > 0)
    Cc <- thiscurve

  mv.av <- function(x,sm) {
    minwindow <- min(nrow(Cc), which(Cc[,1] >= x-(sm/2)))
    maxwindow <- max(1, which(Cc[,1] <= x+(sm/2)))
    mu <- mean(Cc[minwindow:maxwindow,2]) # uniform distribution across the window
    sig <- mean(Cc[minwindow:maxwindow,3])
    return(c(mu,sig))
  }
  av.cc <- sapply(Cc[,1], mv.av, sm=smooth)
  Cc <- cbind(Cc[,1], t(av.cc))

  if(save)
    write.table(Cc, file.path(cc.dir, name), row.names=FALSE, col.names=FALSE, sep=sep)
  invisible(Cc)
}



#' @name calibratable
#' @title Make a table of calibrated dates
#' @description Calibrate a number of radiocarbon dates and make a table containing the calibrated ranges (both highest posterior densities and quantiles).
#' @details Calibration is done taking into account calibration curves and any age offsets (deltaR, deltaSTD). The table will be displayed in an Internet browser, or alternatively saved to a .docx file. Based on ideas by Dr. Vegard Martinsen, NMBU, Norway).
#' If you have a mix of pre- and postbomb dates, they can be calibratabled as follows:
#' y <- c(.9, .8, 1.5); er <- rep(.01, 3)
#' calibrate.table(y, er, is.F=T, thiscurve=glue.ccurves(,1), BCAD=T)
#' @param y The radiocarbon dates
#' @param er The laboratory errors of the radiocarbon dates
#' @param lab The labels of the radiocarbon dates (if any)
#' @param cc The calibration curve to smooth. Calibration curve for 14C dates: 'cc=1' for IntCal20 (northern hemisphere terrestrial), 'cc=2' for Marine20 (marine), 'cc=3' for SHCal20 (southern hemisphere terrestrial). Alternatively, one can also write, e.g., "IntCal20", "Marine13". One can also make a custom-built calibration curve, e.g. using 'mix.ccurves()', and load this using 'cc=4'. In this case, it is recommended to place the custom calibration curve in its own directory, using 'cc.dir' (see below). Explanations of the numbers are provided in the table footer. If there is more than one cc provided, they will be printed in an extra table column.
#' @param BCAD Which calendar scale to use. Defaults to cal BP, \code{BCAD=FALSE}. For the BCAD scale, BC ages are negative.
#' @param postbomb Use 'postbomb=TRUE' to get a postbomb calibration curve (default 'postbomb=FALSE'). For monthly data, type e.g. 'ccurve("sh1-2_monthly")'
#' @param bombalert Warn if a date is close to the lower limit of the IntCal curve. Defaults to \code{postbomb=TRUE}.
#' @param glue Glue postbomb and prebomb curves together. Defaults to 0 (none), can be 1 (IntCal20 + NH1), 2 (IntCal20 + NH2), 3 (IntCal20 + NH3), 4 (SHCal20 + SH1-2) or 5 (SHCal20 + SH3). Note that this will override the value of cc.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., 'cc.dir="ccurves"'.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). Defaults to c().
#' @param is.F Set this to TRUE if the provided age and error are in the F14C timescale.
#' @param is.pMC Set this to TRUE if the provided age and error are in the pMC timescale.
#' @param deltaR Age offset (e.g. for marine samples). If provided, the deltaR and deltaSTD values will be provided as an extra table column.
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation). If provided, the deltaR and deltaSTD values will be provided as an extra table column.
#' @param prob Probability range which should be calculated. Default \code{prob=0.95}.
#' @param age.round Rounding for ages. Defaults to 0 decimals.
#' @param prob.round Rounding for reported probabilities. Defaults to 1 decimal.
#' @param docx By default, the table is written to your web browser. If you wish to write it to a MS-Word document instead, provide the file (with .docx extension) and its location here, e.g., \code{docx="C14_table.docx"}.
#' @author Maarten Blaauw
#' @examples
#'  calibratable(130, 20, bombalert=FALSE)
#'  data(shroud)
#'  calibratable(shroud$y, shroud$er, shroud$ID)
#' @export
calibratable <- function(y, er, lab=c(), cc=1, BCAD=FALSE, postbomb=FALSE, bombalert=TRUE, glue=0, cc.dir=c(), thiscurve=c(), is.F=FALSE, is.pMC=FALSE, deltaR=0, deltaSTD=0, prob=0.95, prob.round=1, age.round=0, docx=c()) {
  if(!requireNamespace("flextable", quietly = TRUE))
    stop("Package 'flextable' is required. Install with install.packages('flextable').")

  if(length(cc) == 1)
    cc <- rep(cc, length(y))
  if(length(er) != length(y))
    stop("er needs to be of the same length as y")
  if(length(cc) != length(y))
    stop("cc needs to be of the same length as y")
  if(length(lab) == 0)
    lab <- 1:length(y) else
      if(length(lab) != length(y))
        stop("lab (labels) needs to be of the same length as y")

  # to be filled later on
  labs <- c(); ys <- c(); ers <- c(); froms <- c(); tos <- c(); percs <- c()
  ccs <- c(); min.rngs <- c(); max.rngs <- c(); grp_idx <- c()

  hasdeltaR <- FALSE
  if(length(deltaR) == 1)
    deltaR <- rep(deltaR, length(y))
  if(length(deltaR) != length(y))
    stop("deltaR needs to be of length 1 or the same length as y")
  if(length(deltaSTD) == 1)
    deltaSTD <- rep(deltaSTD, length(y))
  if(any(deltaR != 0)) {
    hasdeltaR <- TRUE
    deltaRs <- c()
    deltaSTDs <- c()
  }

  for(i in 1:length(y)) {
    this.cal <- caldist(y[i], er[i], cc=cc[i], BCAD=BCAD, postbomb=postbomb, bombalert=bombalert,
      glue=glue, thiscurve=thiscurve, cc.dir=cc.dir, is.F=is.F, is.pMC=is.pMC,
      deltaR=deltaR[i], deltaSTD=deltaSTD[i])
    this.cal <- this.cal[order(this.cal[,1]), ]
    this.hpd <- hpd(this.cal, BCAD=BCAD,
      prob.round=prob.round, age.round=age.round)

    cdf <- cumsum(this.cal[,2] / sum(this.cal[,2])) # ensure it sums to 1
   # p <- this.cal[,2] / sum(this.cal[,2])
   # keep <- p > 0
   # cdf <- cumsum(p[keep])
   # ages <- this.cal[keep, 1]

    qtarget <- c((1-prob)/2, 1-(1-prob)/2) # eat into both edges
    qtarget[1] <- max(qtarget[1], min(cdf))
    qtarget[2] <- min(qtarget[2], max(cdf))
    q <- suppressWarnings(approx(cdf, this.cal[,1], qtarget, rule=2)$y)

   if(BCAD && min(q) > 0) {
     min.rng <- round(q[2], age.round)
     max.rng <- round(q[1], age.round)
   } else {
      min.rng <- round(q[1], age.round)
      max.rng <- round(q[2], age.round)
     }

    n <- length(this.hpd$from)
    froms <- c(froms, as.numeric(this.hpd$from))
    tos <- c(tos, as.numeric(this.hpd$to))
    percs <- c(percs, as.numeric(this.hpd$perc))
    labs <- c(labs, rep(lab[i], n))
    ys <- c(ys, rep(y[i], n))
    ers <- c(ers, rep(er[i], n))
    ccs <- c(ccs, rep(cc[i], n))
    min.rngs <- c(min.rngs, rep(min.rng, n))
    max.rngs <- c(max.rngs, rep(max.rng, n))
    grp_idx <- c(grp_idx, rep(i, n))
    if(hasdeltaR) {
      deltaRs <- c(deltaRs, rep(deltaR[i], n))
      deltaSTDs <- c(deltaSTDs, rep(deltaSTD[i], n))
    }
  }

  # format the values
  C14_sd <- paste(ys, " \u00b1 ", ers)
  hpd_str <- paste0(froms, " \u2013 ", tos, " (", percs, "%)")
  rng_str <- paste0(max.rngs, " \u2013 ", min.rngs)
  
  if(hasdeltaR)
    deltaR_str <- paste(deltaRs, "\u00b1", deltaSTDs)

  collapse_idx <- split(seq_along(grp_idx), grp_idx)
  lab_one <- sapply(collapse_idx, function(idx) labs[idx][1])
  C14_one <- sapply(collapse_idx, function(idx) C14_sd[idx][1])
  cc_one <- sapply(collapse_idx, function(idx) ccs[idx][1])
  rng_one <- sapply(collapse_idx, function(idx) rng_str[idx][1])
  hpd_mult <- sapply(collapse_idx, function(idx) paste(hpd_str[idx], collapse = "\n"))
  if(hasdeltaR)
    deltaR_one <- sapply(collapse_idx, function(idx) deltaR_str[idx][1])
  
  # make table labels
  ylab <- if(is.F) "F14C (\u00b1 sd)" else
    if(is.pMC) "pMC (\u00b1 sd)" else
      "14C BP (\u00b1 sd)"

  if(length(thiscurve) > 0) 
    legend_text <- "Custom calibration curve" else
    if(glue > 0) {
      if(glue %in% 1:3) {
        bomb <- c("NH1", "NH2", "NH3")[glue]
        legend_text <- paste0("Calibration curve: IntCal20 (Reimer et al. 2020) & ",
          bomb, " (Hua et al. 2022)")
      } else if(glue %in% 4:5) {
          bomb <- c("SH1–2", "SH3")[glue - 3]
          legend_text <- paste0("Calibration curve: SHCal20 (Hogg et al. 2020) & ",
            bomb, " (Hua et al. 2022)")
        }
    } else {
      if(postbomb)  {
        cc_labels <- c( 
          "1" = "NH1 (Hua et al. 2022)",
          "2" = "NH2 (Hua et al. 2022)",
          "3" = "NH3 (Hua et al. 2022)",
          "4" = "SH1-2 (Hua et al. 2022)",
          "5" = "SH3 (Hua et al. 2022)")			
	  } else
          cc_labels <- c(
            "1" = "IntCal20 (Reimer et al. 2020)",
            "2" = "Marine20 (Heaton et al. 2020)",
           "3" = "SHCal20 (Hogg et al. 2020)")

      uniq_cc <- unique(cc)
      if(length(uniq_cc) > 1) {
        pairs <- paste0(uniq_cc, " = ", cc_labels[as.character(uniq_cc)])
        legend_text <- paste("cc (calibration curve) legend:",
          paste(pairs, collapse = "; "))
      } else
          legend_text <- paste("Calibration curve:", cc_labels[as.character(uniq_cc)])

      cc_map <- function(v)
        unname(ifelse(as.character(v) %in% names(cc_labels),
          cc_labels[as.character(v)], as.character(v)))
      uniq_cc <- unique(cc)
      if(length(uniq_cc) > 1) {
        pairs <- paste0(uniq_cc, " = ", cc_map(uniq_cc))
        legend_text <- paste0("cc (calibration curve) legend: ", paste(pairs, collapse = "; "))
      } else
        legend_text <- paste("Calibration curve:", cc_map(uniq_cc))
    }
	
  if(BCAD)
    if(min(min.rngs, na.rm=TRUE) < 0)
      legend_text <- paste0(legend_text, ". BC ages are negative")
  hpd.label <- ifelse(BCAD, paste0(100*prob, "% hpd ranges (cal BC/AD)"),
    paste0(100*prob, "% hpd ranges (cal BP)"))
  rng.label <- ifelse(BCAD, paste0(100*prob, "% range (cal BC/AD)"),
    paste0(100*prob, "% range (cal BP)"))

  # Assemble table with one row per date
  astable <- data.frame(labID=lab_one, C14=C14_one, rng=rng_one, stringsAsFactors=FALSE)
  if(hasdeltaR)
    astable <- data.frame(astable, deltaR=deltaR_one, stringsAsFactors=FALSE)
  if(length(unique(cc)) > 1)
    astable <- data.frame(astable, cc=cc_one, stringsAsFactors=FALSE)

  tbl <- data.frame(astable, hpd=hpd_mult, stringsAsFactors=FALSE)
  ft <- flextable::flextable(tbl)
  ft <- flextable::set_header_labels(ft, labID="Lab ID", C14=ylab, rng=rng.label, hpd=hpd.label)
  
  # alignment (only include columns that actually exist)
  cols_to_align <- c("labID", "C14", "rng", "hpd")
  if("deltaR" %in% names(tbl))
    cols_to_align <- c(cols_to_align, "deltaR")
  if("cc" %in% names(tbl))
    cols_to_align <- c(cols_to_align, "cc")
  ft <- flextable::align(ft, j=cols_to_align, align="right", part="body")
  ft <- flextable::align(ft, j=cols_to_align, align="right", part="header")

  ft <- flextable::valign(ft, j=intersect(cols_to_align, c("labID", "C14", "deltaR", "cc", "rng")), valign="top", part="body")
  ft <- flextable::autofit(ft)

  # add footer and return the table
  ft <- flextable::add_footer_lines(ft, values=legend_text)
  
  # save as docx, or to the web browser (default)
  if(length(docx) > 0)
    flextable::save_as_docx(ft, path=file.path(docx)) else
      return(ft)
}
