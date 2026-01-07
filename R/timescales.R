
# from/to	calBP		BCAD		b2k			C14			F14C		pMC			Delta14C
# calBP					calBPtoBCAD	calBPtob2k	calBPtoC14	calBPtoF14C	calBPtopMC	calBPtoDelta14C
# BCAD		BCADtocalP				BCADtob2k	BCADtoC14	BCADtoF14C	BCADtopMC	BCADtoDelta14C
# b2k       b2ktocalBP	b2ktoBCAD				b2ktoC14	b2ktoF14C	b2ktopMC	b2ktoDelta14C
# C14		C14tocalBP	C14toBCAD	C14tob2k				C14toF14C	C14topMC	C14toDelta14C
# F14C		NA			NA			NA			F14CtoC14				F14CtopMC	F14CtoDelta14C
# pMC		NA			NA			NA			pMCtoC14	pMCtoF14C				pMCtoDelta14C
# Delta14C		NA			NA			NA			Delta14CtoC14	Delta14CtoF14C	Delta14CtopMC



#' @name fromto
#' @title translate between timescales
#' @details Upon entering a value and its timescale, this function will find the corresponding values in the other timescales. Note that uncertainties are *not* taken into account, and especially going from C14 BP to cal BP and BC/AD ignores many calibration-related uncertainties. Delta14C values are only reported for entered values on the cal BP or BC/AD scale.
#' @param x The value to be translated into other timescales
#' @param from The timescales of the entered value. Can be "calBP" for cal BP, "BCAD" for BC/AD, "C14" for C14 BP, "F14C" for F14C, or "pMC" for pMC. Delta14C cannot be entered as a value (you could enter the corresponding cal BP or BC/AD ages instead).
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param zero Whether or not zero BC/AD should be included. Defaults to none, FALSE.
#' @param width Width of the righthand plot. Calculated automatically by default (older ages get wider windows).
#' @param digits Rounding of the reported values. Defaults to 0 digits.
#' @param C14.col Colour of the 14C calibration curve. Defaults to semi-transparent blue, \code{C14.col=rgb(0,0,1,.5)}.
#' @param Delta14C.col Colour of the Delta14C curve. Defaults to semi-transparent green, \code{Delta14C.col=rgb(0,.4,0,.4)}.
#' @param ka Whether to use years or ka (thousands of years). Defaults to \code{ka=FALSE}.
#' @param cal.rev Reverse the age axis (right panel). Defaults to TRUE.
#' @param legend.size Size of the font of the legend. Defaults to 0.7 of R's standard size. 
#' @return A plot and output showing the translations into the different timescales.
#' @examples
#'   fromto(0, "BCAD")
#'   fromto(2450, "C14")
#' @export
fromto <- function(x, from="calBP", cc=1, postbomb=1, cc.dir=NULL, thiscurve=NULL, zero=FALSE, width=c(), digits=0, C14.col=rgb(0,0,1,.5), Delta14C.col=rgb(0,.4,0,.4), ka=FALSE, cal.rev=TRUE, legend.size=.7) {
  if(cc==2)
    Cc <- rintcal::ccurve(cc, postbomb=FALSE, cc.dir=cc.dir) else
      if(postbomb)
        Cc <- rintcal::glue.ccurves(cc, postbomb, cc.dir=cc.dir) else
          Cc <- rintcal::ccurve(cc, postbomb, cc.dir=cc.dir)   

  if(grepl("bp", tolower(from))) { # cal BP
    if(x < -69.5 || x > 55e3) # limits postbomb & prebomb
      message("cal BP value lies beyond calibration curve")
    if(x < 0 && postbomb==FALSE)
      message("please provide a postbomb curve")
    calbp <- x
    bcad <- calBPtoBCAD(x, zero)
    c14 <- calBPtoC14(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    f14c <- calBPtoF14C(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    pmc <- calBPtopMC(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    delta14c <- calBPtoDelta14C(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    message(calbp, " cal BP equals ", round(bcad, digits), " cal BC/AD, ", round(c14, digits),
      " 14C BP, ", round(f14c, digits+2), " F14C, ", round(pmc, digits+2),
      " pMC, ", round(delta14c, digits+2), " Delta14C (cc=", cc, ")\n")
  } else
  
  if(grepl("bc", tolower(from))) { # BC/AD
    if(x > 2019.5 || x < -53050)
      message("BC/AD value lies beyond calibration curve")
    if(x > 1950 && postbomb==FALSE)
      message("please provide a postbomb curve")
    bcad <- x
    calbp <- BCADtocalBP(x, zero)
    c14 <- BCADtoC14(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    f14c <- BCADtoF14C(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    pmc <- BCADtopMC(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    delta14c <- BCADtoDelta14C(x, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[1]
    message(bcad, " cal BC/AD equals ", round(calbp, digits), " cal BP, ", round(c14, digits), 
      " 14C BP, ", round(f14c, digits+2), " F14C, ", round(pmc, digits+2),
      " pMC, ", round(delta14c, digits+2), " Delta14C (cc=", cc, ")\n")
  } else
  
  if(grepl("c14", tolower(from))) { # C14
    c14 <- x
    calbp <- unique(C14tocalBP(x))
    bcad <- calBPtoBCAD(calbp)
    f14c <- C14toF14C(x)[1]
    pmc <- C14topMC(x)[1]
    delta14c <- NA # do not calculate since no theta
    message(c14, " 14C BP equals c. ", paste(round(calbp, digits), collapse="/"), " cal BP, c. ", 
      paste(round(bcad, digits), collapse="/"), " cal BC/AD, ",
      round(f14c, digits+2), " F14C, ", round(pmc, digits+2), " pMC\n")
  } else
  
  if(grepl("f", tolower(from))) { # F14C
    f14c <- x
    c14 <- F14CtoC14(x)[1]
    calbp <- unique(C14tocalBP(c14))
    bcad <- calBPtoBCAD(calbp)
    pmc <- 100*x
    delta14c <- NA # do not calculate since no theta
    message(f14c, " F14C equals ", round(pmc, digits+2), " pMC, ",
      round(c14, digits), " 14C BP, c. ",
      paste(round(calbp, digits), collapse="/"), " cal BP, c. ",
      paste(round(bcad, digits), collapse="/"), " cal BC/AD\n")
  } else
  
  if(grepl("p", tolower(from))) { # pMC
    pmc <- x
    f14c <- pmc/100
    c14 <- pMCtoC14(x)[1]
    calbp <- unique(C14tocalBP(c14))
    bcad <- calBPtoBCAD(calbp)
    pmc <- 100*x
    delta14c <- NA # do not calculate since no theta
    message(pmc, " pMC equals ", round(f14c, digits), " F14C, ",
      round(c14, digits), " 14C BP, c. ",
      paste(round(calbp, digits), collapse="/"), " cal BP, c. ",
      paste(round(bcad, digits), collapse="/"), " cal BC/AD\n")
  } else
      stop("Please provide a correct entry for parameter 'from'")
  
  if(length(width) == 0) {
    exp_interp <- function(x) {
      b <- log(10e3 / 100) / 55000
      return(100 * exp(b * x))
    }
    width <- max(exp_interp(calbp))
  }

  layout(matrix(1:2, nrow=1), widths=c(.5, .5))
  c14lim <- c(c14-width, c14+width) # hmmm
  c14lim <- range(Cc[,2])
  flim <- c(0, min(1.2*f14c,2))
  par(mar=c(4,3.2,3,1), mgp=c(2, .7, 0), yaxt="s")

  calbp <- as.numeric(calbp)
  f14c <- as.numeric(f14c)
  c14 <- as.numeric(c14)
  plot(f14c, c14, xlim=flim, ylim=c14lim, 
    xlab=expression("F"^14*"C"), ylab=expression(""^14*"C BP"),
      bty="c", type="n", xaxs="i")
  abline(h=c14, lty=2)
  abline(v=f14c, lty=2)
  curve(-8033*log(x), from=min(flim), n=1e3, to=1.05*max(flim), add=TRUE, col=C14.col, lwd=2)
  pr <- pretty(seq(par("usr")[1], par("usr")[2], length=5))
  axis(3, pr, 100*pr, line=0)
  mtext("pMC", 3, 2)

  values <- c(bcad, calbp, c14, f14c, pmc, delta14c)
  names(values) <- c("BC/AD", "cal BP", "C14", "F14C", "pMC", "Delta14C")

  if(length(calbp) == 1)
    calstring <- paste0(round(calbp, digits), " cal BP\n", round(bcad, digits), " cal BC/AD") else
      calstring <- paste0(paste0(round(calbp, digits), collapse="/"), " cal BP\n",
        paste0(round(bcad, digits), collapse="/"), " cal BC/AD")
  c14string <- bquote(.(paste0(round(c14, digits))) ~ {}^{14}*C)
  f14cstring <- bquote(.(paste0(round(f14c, digits+2))) ~ F^{14}*C)
  pmcstring <- paste(round(pmc, digits+1), "pMC")
  delta14cstring <- bquote(.(paste0(round(delta14c, digits+2))) ~ Delta^{14}*C)

  legend.coors <- legend("topleft", legend=1:6, ncol=1, cex=legend.size, xjust=1, plot=FALSE)$text # legend coordinates
  text(legend.coors$x[1], legend.coors$y[1], calstring, col=1, cex=legend.size, adj=0)
  text(legend.coors$x[1], legend.coors$y[3], c14string, col=C14.col, cex=legend.size, adj=0)
  text(legend.coors$x[1], legend.coors$y[4], f14cstring, col=C14.col, cex=legend.size, adj=0)
  text(legend.coors$x[1], legend.coors$y[5], pmcstring, col=C14.col, cex=legend.size, adj=0)
  if(!is.na(delta14c))
    text(legend.coors$x[1], legend.coors$y[6], delta14cstring, col=Delta14C.col, cex=legend.size, adj=0)

  par(mar=c(4,3,3,3))
  mincalbp <- min(calbp) - width
  maxcalbp <- max(calbp) + width
  Delta14C.coors <- draw.ccurve(mincalbp, maxcalbp, cc1=cc, cc2=cc, cc1.postbomb=TRUE, cc2.postbomb=TRUE, timescale2="d", cc2.col=Delta14C.col, cc2.fill=Delta14C.col, add.yaxis=TRUE, cc.dir=cc.dir, ka=ka, bty="n", legend=NA, xaxs="i", yaxt="n", c14.lab="", cal.rev=cal.rev)
  C14.coors <- par("usr")
  if(!is.na(delta14c))
    f.y <- as.numeric(((delta14c-Delta14C.coors[3])/(Delta14C.coors[4] - Delta14C.coors[3])) * (C14.coors[4] - C14.coors[3]) + C14.coors[3])

  pr <- pretty(seq(C14.coors[1], C14.coors[2], length=5))
  if(ka) {
    mtext("k cal BC/AD", 3, 2)
    calbp <- calbp/1e3
    c14 <- c14/1e3
    axis(3, pr, calBPtoBCAD(1e3*pr)/1e3)
  } else {
      mtext("cal BC/AD", 3, 2)
      axis(3, pr, calBPtoBCAD(pr))
    }
  segments(C14.coors[1], C14.coors[3], C14.coors[2], C14.coors[3])
  axis(2, pretty(C14.coors[3:4]), col=C14.col, col.axis=C14.col)
  mtext(expression(""^14*"C BP"), 2, 2, col=C14.col)
  abline(v=C14.coors[2], col=Delta14C.col)
  mtext(expression(Delta^14*C), 4, 2, col=Delta14C.col)

  abline(v=calbp, lty=2)
  segments(C14.coors[1]-1e3, c14, calbp, c14, lty=2, col=C14.col)
  if(!is.na(delta14c))
    segments(calbp, f.y, C14.coors[2]+1e3, f.y, lty=2, col=Delta14C.col)

  invisible(values)
}



#' @name calBPtoBCAD
#' @title calculate cal BC/AD ages from cal BP ages
#' @details Turn cal BP ages into cal BC/AD (equivalent to cal BCE/CE). Negative ages indicate cal BC, positive ages cal AD. Since the Gregorian and Julian calendars do not include 0 BCAD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from -1 (i.e., 1 BC) to 1 AD. Other calendars, such as the astronomical one, do include zero. 
#' @param x The calBP age(s) to be translated into cal BC/AD ages. 
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @return The cal BC/AD age(s). BC ages are negative, AD ages are positive.
#' @examples
#'  calBPtoBCAD(2024)
#'  calBPtoBCAD(1945:1955, zero=FALSE)
#'  calBPtoBCAD(1945:1955, zero=TRUE)
#' @export 
calBPtoBCAD <- function(x, zero=FALSE) {
  xv <- as.numeric(x)

  if(zero) 
    X <- 1950 - xv else {
      X <- 1950 - xv
      idx <- xv >= 1950
      if(any(idx)) {
        X[idx] <- 1949 - xv[idx]
      }
    }

  dim(X) <- dim(x)
  dimnames(X) <- dimnames(x)

  return(X)  
}



#' @name calBPtob2k
#' @title calculate b2k ages from cal BP ages
#' @details Turn cal BP ages into b2k ages (years before AD 2000), which are often used in the ice core community.
#' @param x The calBP age(s) to be translated into b2k ages.
#' @return The b2k ages.
#' @examples
#'  calBPtob2k(-50)
#' @export
calBPtob2k <- function(x)
  return(x + 50)



#' @name calBPtoC14
#' @title Find the 14C age and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding 14C age and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param x The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#' calBPtoC14(100)
#' @export
calBPtoC14 <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  if(is.data.frame(x) && ncol(x) == 1)
    x <- x[[1]] # accept either a data.frame or a vector
  if(is.null(thiscurve)) {
    if(cc == 2) # Marine20 has no postbomb counterpart
      cc <- rintcal::ccurve(cc=cc, postbomb=postbomb, cc.dir=cc.dir) else
        if(postbomb)
          cc <- rintcal::glue.ccurves(prebomb=cc, postbomb=postbomb, cc.dir=cc.dir) else
            cc <- rintcal::ccurve(cc=cc, postbomb=postbomb, cc.dir=cc.dir) 
    } else
        cc <- thiscurve
  mu <- approx(cc[,1], cc[,2], x, rule=rule)$y
  er <- approx(cc[,1], cc[,3], x, rule=rule)$y

  return(round(data.frame(C14=mu, er=er), roundby))
}



#' @name calBPtoF14C
#' @title Find the F14C and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param x The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   calBPtoF14C(100)
#' @export
calBPtoF14C <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  y <- calBPtoC14(x, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=Inf)
  return(C14toF14C(y[,1], y[,2], roundby=roundby)) # is a data.frame; any rounding done at the end
}



#' @name calBPtopMC
#' @title Find the pMC and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param x The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   calBPtopMC(100)
#' @export
calBPtopMC <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  y <- calBPtoC14(x, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=Inf)
  return(C14topMC(y[,1], y[,2], roundby=roundby)) # data.frame; any rounding done at the end
}



#' @name calBPtoDelta14C
#' @title Find the Delta14C and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding Delta14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve Delta14C belonging to the entered cal BP age
#' @param x The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   calBPtoDelta14C(100)
#' @export
calBPtoDelta14C <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  F <- calBPtoF14C(x, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=Inf)
  Deltamn <- F14CtoDelta14C(F[,1], t=x)
  Deltaup <- F14CtoDelta14C(F[,1]+F[,2], t=x)
  sdev <- Deltaup-Deltamn

  return(round(data.frame(Delta14C=Deltamn, er=sdev), roundby))
}



#' @name BCADtocalBP
#' @title calculate cal BP ages from cal BC/AD ages
#' @details Turn cal BC/AD (or BCE/CE) ages into cal BP ages. Negative ages indicate BC, positive ages AD. Since the Gregorian and Julian calendars do not include 0 BC/AD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from -1 (i.e., 1 BC) to 1 AD. Other calendars, such as the astronomical one, do include zero. The often-used BCE/CE ages are equivalent to BC/AD.
#' @param x The cal BCAD age(s) to be translated into cal BP age(s). BC ages are negative, AD ages are positive.
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @return The cal BP age(s).
#' @examples
#'  BCADtocalBP(2025)
#'  BCADtocalBP(-1, zero=TRUE)
#'  BCADtocalBP(-1, zero=FALSE)
#' @export 
BCADtocalBP <- function(x, zero=FALSE) 
  if(zero)
    return(1950-x) else {
      X <- 1950-x 
      neg <- which(x < 0)  
      if(length(neg) > 0)
        X[neg] <- 1949-x[neg]
        return(X) # X = 0 is not an error
    }



#' @name BCADtob2k
#' @title calculate b2k from cal BC/AD ages
#' @details Turn cal BC/AD (or BCE/CE) ages into b2k ages. b2k ages are used frequently in the ice core community. Negative ages indicate BC, positive ages AD. Since the Gregorian and Julian calendars do not include 0 BC/AD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from -1 (i.e., 1 BC) to 1 AD. Other calendars, such as the astronomical one, do include zero. The often-used BCE/CE ages are equivalent to BC/AD.
#' @param x The BCAD age(s) to be translated into b2k age(s). BC ages are negative, AD ages are positive.
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @return The b2k age(s).
#' @examples
#'  BCADtob2k(2025)
#'  BCADtob2k(-1, zero=TRUE)
#'  BCADtob2k(-1, zero=FALSE)
#' @export
BCADtob2k <- function(x, zero=FALSE)
  if(zero)
    return(2000-x) else {
      X <- 2000-x
      neg <- which(x < 0)
      if(length(neg) > 0)
        X[neg] <- 1999-x[neg]
        return(X) # x = 0 is not an error
    }



#' @name BCADtoC14
#' @title Find the 14C age and error belonging to a cal BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding 14C age and error are returned. BC ages are negative. In this implementation, the year 0 BC/AD does exist.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than AD 1950, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered BC/AD age
#' @param x The cal BC/AD year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   BCADtoC14(100)
#' @export
BCADtoC14 <- function(x, cc=1, postbomb=FALSE, zero=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf)
  return(calBPtoC14(BCADtocalBP(x, zero=zero),
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=roundby))



#' @name BCADtoF14C
#' @title Find the F14C and error belonging to a cal BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C and error are returned. BC ages are negative.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than AD 1950, a postbomb curve will have to be provided.
#' @return The calibration-curve F14C belonging to the entered BC/AD age
#' @param x The cal BC/AD year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   BCADtoF14C(100)
#' @export
BCADtoF14C <- function(x, cc=1, postbomb=FALSE, zero=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf)
  return(calBPtoF14C(BCADtocalBP(x, zero=zero), 
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=roundby))



#' @name BCADtopMC
#' @title Find the pMC and error belonging to a cal BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding pMC and error are returned. BC ages are negative.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than AD 1950, a postbomb curve will have to be provided.
#' @return The calibration-curve F14C belonging to the entered cal BC/AD age
#' @param x The cal BC/AD year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   BCADtopMC(100)
#' @export
BCADtopMC <- function(x, cc=1, postbomb=FALSE, zero=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf)
  return(calBPtopMC(BCADtocalBP(x, zero=zero), 
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=roundby))



#' @name BCADtoDelta14C
#' @title Find the Delta14C and error belonging to a cal BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding Delta14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BC/AD age.
#' @param x The cal BC/AD year.
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @param cc calibration curve (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   BCADtoDelta14C(1900)
#' @export
BCADtoDelta14C <- function(x, zero=FALSE, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  calBP <- BCADtocalBP(x, zero)
  Fres <- calBPtoF14C(calBP, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=Inf)
  Deltamn <- F14CtoDelta14C(Fres[,1], t=calBP)
  Deltaup <- F14CtoDelta14C(Fres[,1]+Fres[,2], t=calBP)

  return(round(data.frame(Delta14C=Deltamn, er=Deltaup-Deltamn), roundby))
}



#' @name b2ktocalBP
#' @title calculate cal BP ages from b2k ages
#' @details Turn b2k ages (often used in the ice core community, AD 2000) into cal BP ages.
#' @param x The b2k age(s) to be translated into cal BP age(s).
#' @return The cal BP age(s).
#' @examples
#'  b2ktocalBP(0)
#' @export
b2ktocalBP <- function(x)
  return(x - 50)



#' @name b2ktoBCAD
#' @title calculate cal BC/AD ages from b2k ages
#' @details Turn b2k ages (popular in the ice core community) into cal BC/AD (or cal BCE/CE). Negative ages indicate BC, positive ages AD. Since the Gregorian and Julian calendars do not include 0 BCAD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from -1 (i.e., 1 BC) to 1 AD. Other calendars, such as the astronomical one, do include zero. The often-used BCE/CE ages are equivalent to BC/AD.
#' @param x The b2k age(s) to be translated into cal BC/AD ages.
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @return The cal BC/AD age(s). BC ages are negative, AD ages are positive.
#' @examples
#'  b2ktoBCAD(0)
#'  b2ktoBCAD(1990:2010, zero=TRUE)
#'  b2ktoBCAD(1990:2010, zero=FALSE)
#' @export
b2ktoBCAD <- function(x, zero=FALSE) {
  if(zero)
    return(2000-x) else {
      X <- 2000-x
      neg <- which(x >= 2000)
      if(length(neg) > 0)
        X[neg] <- 1999 - x[neg]
      return(X) # x = 0 is not an error
    }
}



#' @name b2ktoC14
#' @title Find the 14C age and error belonging to a b2k age.
#' @description Given a b2k age (years before AD 2000, popular in the ice core community), the calibration curve (default cc=1) is interpolated and the corresponding 14C age and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than 50 b2k, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered b2k age
#' @param x The b2k year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   b2ktoC14(100)
#' @export
b2ktoC14 <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf)
  return(calBPtoC14(b2ktocalBP(x),
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=roundby))



#' @name b2ktoF14C
#' @title Find the F14C and error belonging to a b2k age.
#' @description Given a b2k age (years before AD 2000, popular in the ice core community), the calibration curve (default cc=1) is interpolated and the corresponding F14C and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than 50 b2k, a postbomb curve will have to be provided.
#' @return The calibration-curve F14C belonging to the entered b2k age
#' @param x The b2k year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   b2ktoF14C(100)
#' @export
b2ktoF14C <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf)
  return(calBPtoF14C(b2ktocalBP(x),
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=roundby))



#' @name b2ktopMC
#' @title Find the pMC and error belonging to a b2k age.
#' @description Given a b2k age (years before AD 2000, popular in the ice core community), the calibration curve (default cc=1) is interpolated and the corresponding pMC and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than 50 b2k, a postbomb curve will have to be provided.
#' @return The calibration-curve F14C belonging to the entered b2k age
#' @param x The b2k year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   b2ktopMC(100)
#' @export
b2ktopMC <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf)
  return(calBPtopMC(b2ktocalBP(x),
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=roundby))



#' @name b2ktoDelta14C
#' @title Find the Delta14C and error belonging to a b2k age.
#' @description Given a b2k age (years before AD 2000, popular in the ice core community), the calibration curve (default cc=1) is interpolated and the corresponding Delta14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For b2k < 50, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered b2k age
#' @param x The b2k year.
#' @param cc calibration curve (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   b2ktoDelta14C(100)
#' @export
b2ktoDelta14C <- function(x, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  calBP <- b2ktocalBP(x)
  Fres <- calBPtoF14C(calBP, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve, roundby=Inf)
  Deltamn <- F14CtoDelta14C(Fres[,1], t=calBP)[,1]
  Deltaup <- F14CtoDelta14C(Fres[,1]+Fres[,2], t=calBP)[,1]

  return(round(data.frame(Delta14C=Deltamn, er=Deltaup-Deltamn), roundby))
}



#' @name C14tocalBP
#' @title Find the calBP age(s) crossing a C14 age.
#' @description Find the cal BP ages where the calibration curve crosses a given C14 age. This function is for illustration only and not to be used for, e.g., calibration, because intercept calibration is an outdated method.
#' @details. Whereas each cal BP age will only have one single IntCal radiocarbon age (mu), the same cannot be said for the other way round. Recurring C14 ages do happen, especially during periods of plateaux and wiggles. Therefore, there can be multiple cal BP ages for a single C14 age. In the early days, radiocarbon calibration used an 'intercept method' to find possible calendar ages belonging to a radiocarbon age, but this is problematic since small deviations in the C14 age can easily cause more or fewer crossing cal BP ages (try for example C14tocalBP(130) vs C14tocalBP(129)), and moreover, this approach does not deal well with the errors in either a date of the calibration curve. Therefore, the probabilistic approach to radiocarbon calibration (which starts with a cal BP age and then looks up the corresponding C14 age) has taken over as the standard. 
#' @return The cal BP age(s) belonging to the entered C14 age
#' @param y The C14 age. No errors are assumed. Can only deal with one C14 age at a time.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   y <- 130
#'   calibrate(y,10)
#'   abline(h=y)
#'   abline(v=C14tocalBP(y))
#' @export
C14tocalBP <- function(y, cc=1, postbomb=FALSE, rule=2, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  if(length(y) != 1)
    stop("can only deal with one y at a time")
  if(is.null(thiscurve)) {
    if(cc == 2) # Marine20 has no postbomb counterpart
      cc <- rintcal::ccurve(cc=cc, postbomb=postbomb, cc.dir=cc.dir) else
        if(postbomb)
          cc <- rintcal::glue.ccurves(prebomb=cc, postbomb=postbomb, cc.dir=cc.dir) else
            cc <- rintcal::ccurve(cc=cc, postbomb=postbomb, cc.dir=cc.dir)
    } else
        cc <- thiscurve

  sel <- which((cc[,2] < (y+100) ) * (cc[,2] > (y-100)) > 0)
  if(length(sel) == 0) # y lies beyond the curve
    return(NA)

  cc <- cc[sel,]
  x <- seq(min(cc[,1]), max(cc[,1]), length=1e3)
  mu <- approx(cc[,1], cc[,2], x)$y

  # Identify indices where y crosses the threshold
  sign_changes <- which((mu[-1] - y) * (mu[-length(mu)] - y) < 0)

  # Perform linear interpolation to estimate exact crossing points
  intercepts <- c()
  for(i in seq_along(sign_changes)) {
    idx <- sign_changes[i]
    x1 <- x[idx]; x2 <- x[idx + 1]
    mu1 <- mu[idx]; mu2 <- mu[idx + 1]
    intercepts[i] <- x1 + (x2 - x1) * (y - mu1) / (mu2 - mu1)
  }
  return(round(intercepts, roundby))
}



#' @name C14toBCAD
#' @title Find the cal BCAD age(s) crossing a C14 age.
#' @description Find the BCAD ages where the calibration curve crosses a given C14 age. This function is for illustration only and not to be used for, e.g., calibration, because intercept calibration is an outdated method.
#' @return The cal BCAD age(s) belonging to the entered C14 age
#' @details. Whereas each cal BC/AD age will only have one single IntCal radiocarbon age (mu), the same cannot be said for the other way round. Recurring C14 ages do happen, especially during periods of plateaux and wiggles. Therefore, there can be multiple cal BC/AD ages for a single C14 age. In the early days, radiocarbon calibration used an 'intercept method' to find possible calendar ages belonging to a radiocarbon age, but this is problematic since small deviations in the C14 age can easily cause more or fewer crossing cal BC/AD ages (try for example C14toBCAD(130) vs C14toBCAD(129)), and moreover, this approach does not deal well with the errors in either a date of the calibration curve. Therefore, the probabilistic approach to radiocarbon calibration (which starts with a cal BC/AD age and then looks up the corresponding C14 age) has taken over as the standard.
#' @param y The C14 age.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param zero Whether or not zero BC/AD should be included. Defaults to \code{zero=FALSE}.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'   y <- 130
#'   calibrate(y, 10, BCAD=TRUE)
#'   abline(h=y)
#'   abline(v=C14toBCAD(y))
#' @export
C14toBCAD <- function(y, cc=1, postbomb=FALSE, rule=1, zero=FALSE, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  x <- C14tocalBP(y, roundby=Inf)
  return(round(calBPtoBCAD(x, zero=zero), roundby))
}



#' @name C14tob2k
#' @title Find the b2k age(s) crossing a C14 age.
#' @description Find the b2k ages (years before AD 2000, popular in the ice core community) where the calibration curve crosses a given C14 age. This function is for illustration only and not to be used for, e.g., calibration, because intercept calibration is an outdated method.
#' @return The b2k age(s) belonging to the entered C14 age
#' @details. Whereas each calendar age will only have one single IntCal radiocarbon age (mu), the same cannot be said for the other way round. Recurring C14 ages do happen, especially during periods of plateaux and wiggles. Therefore, there can be multiple cal BP ages for a single C14 age. In the early days, radiocarbon calibration used an 'intercept method' to find possible calendar ages belonging to a radiocarbon age, but this is problematic since small deviations in the C14 age can easily cause more or fewer crossing cal BP ages (try for example C14tocalBP(130) vs C14tocalBP(129)), and moreover, this approach does not deal well with the errors in either a date of the calibration curve. Therefore, the probabilistic approach to radiocarbon calibration (which starts with a cal BP age and then looks up the corresponding C14 age) has taken over as the standard.
#' @param y The C14 age.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @author Maarten Blaauw
#' @examples
#'  C14tob2k(130, 20)
#' @export
C14tob2k <- function(y, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL, roundby=Inf) {
  x <- C14tocalBP(y, roundby=Inf)
  return(round(calBPtob2k(x),roundby))
}



#' @name C14toF14C
#' @title Calculate F14C values from C14 ages
#' @description Calculate F14C values from radiocarbon ages
#' @details Post-bomb dates are often reported as F14C or fraction modern carbon. Since software such as Bacon expects radiocarbon ages,
#' this function can be used to calculate F14C values from radiocarbon ages. The reverse function of \link{F14C.age}.
#' @param y Reported mean of the 14C age.
#' @param er Reported error of the 14C age. If left empty, will translate y to F14C.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years).
#' @param botherrors Since going from C14 to F14C involves a logarithmic transformation (F=exp(-y/lambda)), errors that are symmetric on the C14 scale will become asymmetric on the F14C scale. By default, only the largest error is reported, but if \code{botherrors=TRUE}, both errors are reported.
#' @return F14C values from C14 ages.
#' @examples
#'   C14toF14C(-2000, 20)
#' @export
C14toF14C <- function(y, er=NULL, roundby=Inf, lambda=8033, botherrors=FALSE) {
  y <- as.matrix(y)

  if(!is.null(er)) {
    if(is.vector(er))
      er <- matrix(er, ncol = 1)
    if(is.vector(y))
      y <- matrix(y, ncol = 1)
    if(nrow(er) != nrow(y))
      stop("y and er must have same length")
  }

  if(ncol(y) != 1 || (!is.null(er) && ncol(er) != 1))
    stop("y and er must be single-column vectors or matrices")

  fy <- exp(-y / lambda)

  if(is.null(er))
    return(round(fy, roundby))

  er1 <- abs(fy - exp(-(y - er) / lambda)) # younger error
  er2 <- abs(fy - exp(-(y + er) / lambda)) # older error

  if(botherrors) {
    result <- data.frame(F14C=as.vector(fy), er1=as.vector(er1), er2=as.vector(er2))
  } else {
      result <- data.frame(F14C=as.vector(fy), er=as.vector(pmax(er1, er2)))
    }
  return(round(result, roundby))
}



#' @name CtoF
#' @title Calculate F14C values from C14 ages
#' @description Calculate F14C values from radiocarbon ages
#' @details Post-bomb dates are often reported as F14C or fraction modern carbon. Since software such as Bacon expects radiocarbon ages,
#' this function can be used to calculate F14C values from radiocarbon ages. The reverse function of \link{F14C.age}. This function is a shortcut to C14toF14C.
#' @param y Reported mean of the 14C age.
#' @param er Reported error of the 14C age. If left empty, will translate y to F14C.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years).
#' @param botherrors Since going from C14 to F14C involves a logarithmic transformation (F=exp(-y/lambda)), errors that are symmetric on the C14 scale will become asymmetric on the F14C scale. By default, only the largest error is reported, but if \code{botherrors=TRUE}, both errors are reported.
#' @return F14C values from C14 ages.
#' @examples
#'   CtoF(-2000, 20)
#' @export
CtoF <- function(y, er=NULL, roundby=Inf, lambda=8033, botherrors=FALSE)
  C14toF14C(y=y, er=er, roundby=roundby, lambda=lambda, botherrors=botherrors)



#' @name C14topMC
#' @title Calculate pMC values from C14 ages
#' @description Calculate pMC values from radiocarbon ages
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#' this function can be used to calculate pMC values from radiocarbon ages. The reverse function of \link{pMCtoC14}.
#' @param y Reported mean of the C14 age.
#' @param er Reported error of the C14 age.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return pMC values from C14 ages.
#' @examples
#'   C14topMC(-2000, 20)
#'   C14topMC(-2000, 20, 1)
#' @export
C14topMC <- function(y, er=NULL, roundby=Inf, lambda=8033) {
  return(round(100*C14toF14C(y=y, er=er, roundby=Inf, lambda=lambda), roundby))
#  if(length(er)==0)
#    return(round(data.frame(pMC=100*toF[,1]), roundby)) else
#      return(round(data.frame(pMC=100*toF[,1], sd=100*toF[,2]), roundby))
}



#' @name C14toDelta14C
#' @title Transform C14 age(s) into Delta14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates C14 ages into Delta14C, the historical level of Delta14C in the year t cal BP. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param y The C14 age to translate
#' @param er Reported error of the C14 age. Returns just the mean if left empty.
#' @param t the cal BP age
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda Radiocarbon's mean-life, based on the Cambridge half-life
#' @return The corresponding Delta14C value
#' @examples
#'   C14toDelta14C(0.985, 20, 222)
#' @export
C14toDelta14C <- function(y, er=NULL, t, roundby=Inf, lambda=5730/log(2)) {
  if(!length(y) == length(t))
    stop("inputs 'y' and 't' must have the same length")
  asF <- cbind(C14toF14C(cbind(y), cbind(er), roundby=Inf))
  Deltamn <- 1000 * ((asF[,1] / exp(cbind(-1*t)/lambda)) - 1)
  if(is.null(er)) {
    return(round(data.frame(Delta14C=Deltamn), roundby))
  } else {
      Fup <- asF[,1] + asF[,2]
      Fdown <- asF[,1] - asF[,2]
      Deltaup <- 1000 * ((Fup / exp(-t / lambda)) - 1)
      Deltadown <- 1000 * ((Fdown / exp(-t / lambda)) - 1)
      sdev <- pmax(abs(Deltaup - Deltamn), abs(Deltadown - Deltamn))
      return(round(data.frame(Delta14C=Deltamn, er=sdev), roundby))
    }
}



#' @name F14CtoC14
#' @title Calculate C14 ages from F14C values.
#' @description Calculate C14 ages from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param F14C Reported mean of the F14C
#' @param er Reported error of the F14C. Returns just the mean if left empty.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years).
#' @param botherrors Since going from C14 to F14C involves a logarithmic transformation (F=exp(-y/lambda)), errors that are symmetric on the C14 scale will become asymmetric on the F14C scale. By default, only the largest error is reported, but if \code{reportbotherrors=TRUE}, both errors are reported. #' @return The radiocarbon ages from the F14C values. If F14C values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   F14CtoC14(1.10, 0.005) # a postbomb date, so with a negative C14 age
#'   F14CtoC14(.80, 0.005) # prebomb dates can also be calculated
#' @export
F14CtoC14 <- function(F14C, er=NULL, roundby=Inf, lambda=8033, botherrors=FALSE) {
  y <- rep(NaN, length(F14C))
  valid.y <- which(F14C > 0)
  y[valid.y] <- -lambda * log(F14C[valid.y])
  y[F14C == 0] <- Inf

  if(is.null(er)) # no errors
    return(round(y, roundby))
  
  if(min(er) < 0)
    stop("cannot have negative errors")
  
  # report the older (=larger) error if only reporting one
  valid.er <- which(F14C - er > 0)
  error.older <- rep(NaN, length(y))
  error.older[valid.er] <-
    abs(y[valid.er] - (-lambda * log(F14C[valid.er] - er[valid.er])))

  if(botherrors) {
    error.younger <- rep(NaN, length(y))
    error.younger[valid.y] <- y[valid.y] - (-lambda * log(F14C[valid.y] + er[valid.y]))
    result <- data.frame(C14=y, er.younger=error.younger, er.older=error.older)
  } else
      result <- data.frame(C14=y, er=error.older)

  return(round(result, roundby))
}



#' @name FtoC
#' @title Calculate C14 ages from F14C values.
#' @description Calculate C14 ages from F14C values of radiocarbon dates. Shorthand for the function F14CtoC14.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param F14C Reported mean of the F14C
#' @param er Reported error of the F14C. Returns just the mean if left empty.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years).
#' @param botherrors Since going from C14 to F14C involves a logarithmic transformation (F=exp(-y/lambda)), errors that are symmetric on the C14 scale will become asymmetric on the F14C scale. By default, only the largest error is reported, but if \code{reportbotherrors=TRUE}, both errors are reported. #' @return The radiocarbon ages from the F14C values. If F14C values are above 100\%, the resulting radiocarbon ages will be negative.
#' @return The radiocarbon ages from the F14C values. If F14C values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   FtoC(1.10, 0.005) # a postbomb date, so with a negative C14 age
#'   FtoC(.80, 0.005) # prebomb dates can also be calculated
#' @export
FtoC <- function(F14C, er=NULL, roundby=Inf, lambda=8033, botherrors=FALSE) 
  F14CtoC14(F14C=F14C, er=er, roundby=roundby, lambda=lambda, botherrors=botherrors)



#' @name F14CtopMC
#' @title Calculate pMC ages from F14C values.
#' @description Calculate pMC values from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param F14C Reported mean of the F14C
#' @param er Reported error of the F14C. Returns just the mean if left empty.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @return The pMC values from the F14C values. Basically the original values multiplied by 100.
#' @examples
#'   F14CtopMC(1.10, 0.5)
#' @export
F14CtopMC <- function(F14C, er=NULL, roundby=Inf)
  return(round(data.frame(pMC=100*F14C, er=100*er), roundby))



#' @name F14CtoDelta14C
#' @title Transform F14C into Delta14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates F14C values into Delta14C, the historical level of Delta14C in the year t cal BP. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param F14C The F14C value to translate
#' @param er Reported error of the F14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda Radiocarbon's mean-life, based on the Cambridge half-life
#' @return The corresponding Delta14C value
#' @examples
#'   F14CtoDelta14C(0.89, .001, 900)
#' @export
F14CtoDelta14C <- function(F14C, er=NULL, t, roundby=Inf, lambda=5730/log(2)) {
  # If er is NULL, return the Deltamn for F14C and t
  Deltamn <- 1000 * ((F14C / exp(-t / lambda)) - 1)
  if(is.null(er)) {
    return(round(data.frame(Delta14C=Deltamn), roundby))
  } else {
      Deltaup <- 1000 * (((F14C + er) / exp(-t / lambda)) - 1)
      return(round(data.frame(Delta14C=Deltamn, er=Deltaup - Deltamn), roundby))
    }
}



#' @name pMCtoC14
#' @title Calculate C14 ages from pMC values.
#' @description Calculate C14 ages from pMC values of radiocarbon dates.
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from pMC values. The reverse function is C14.pMC.
#' @param pMC Reported mean of the pMC.
#' @param er Reported error of the pMC.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return Radiocarbon ages from pMC values. If pMC values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   pMCtoC14(110, 0.5) # a postbomb date, so with a negative 14C age
#'   pMCtoC14(80, 0.5) # prebomb dates can also be calculated
#'   pMCtoC14(.8, 0.005) # throws a warning, use F14C.age instead
#' @export
pMCtoC14 <- function(pMC, er=NULL, roundby=Inf, lambda=8033) {
  y <- -lambda * log(pMC/100)
  if(is.null(er)) 
    return(round(y, roundby))
    else {
      sdev <- y - -lambda * log((pMC+er)/100)
      return(round(data.frame(C14=y, er=sdev), roundby))
    }
}



#' @name pMCtoF14C
#' @title Calculate pMC ages from F14C values.
#' @description Calculate pMC values from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param pMC Reported mean of the F14C
#' @param er Reported error of the pMC value. Returns just the mean if left empty.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @return The F14C values from the pMC values. Basically the original values divided by 100.
#' @examples
#'   pMCtoF14C(110, 5)
#' @export
pMCtoF14C <- function(pMC, er=NULL, roundby=Inf)
  if(length(er) == 0)
    return(round(data.frame(F14C=pMC/100), roundby)) else 
      return(round(data.frame(F14C=pMC/100, er=er/100), roundby))



#' @name pMCtoDelta14C
#' @title Transform pMC into Delta14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates F14C values into Delta14C, the historical level of Delta14C in the year t cal BP. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param pMC The pMC value to translate
#' @param er Reported error of the pMC value. Returns just the mean if left empty.
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param t the cal BP age
#' @return The corresponding Delta14C value
#' @examples
#'   pMCtoDelta14C(98.5, 1, 222)
#' @export
pMCtoDelta14C <- function(pMC, er=NULL, t, roundby=Inf) {
  asF <- pMCtoF14C(pMC, er, roundby=Inf)
  if(length(er) == 0)
    return(round(F14CtoDelta14C(asF[,1], NULL, t), roundby)) else
      return(round(F14CtoDelta14C(asF[,1], asF[,2], t), roundby)) 
}



#' @name Delta14CtoC14
#' @title Transform Delta14C into C14 age
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates Delta14C, the historical level of Delta14C in the year t cal BP, to C14 ages. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param Delta14C The Delta14C value to translate
#' @param er Reported error of the Delta14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @return The corresponding C14 age
#' @examples
#'   Delta14CtoC14(-10, 1, 238)
#' @export
Delta14CtoC14 <- function(Delta14C, er=NULL, t, roundby=Inf) {
  toF <- Delta14CtoF14C(Delta14C=Delta14C, er=er, t=t, roundby=Inf)
  if(ncol(toF) == 1) {
    return(C14=F14CtoC14(toF[,1], c(), roundby=roundby))
  } else {
      toF <- as.matrix(toF, ncol=2)
      return(F14CtoC14(toF[,1], toF[,2], roundby=roundby))
    }
}



#' @name Delta14CtoF14C
#' @title Transform Delta14C into F14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates Delta14C, the historical level of Delta14C in the year t cal BP, to F14C values. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param Delta14C The Delta14C value to translate
#' @param er Reported error of the Delta14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @param lambda Radiocarbon's mean-life, based on the Cambridge half-life
#' @return The corresponding F14C value
#' @examples
#'   Delta14CtoF14C(-10, 1, 238)
#' @export
Delta14CtoF14C <- function(Delta14C, er=NULL, t, roundby=Inf, lambda=5730/log(2)) {
  asF <- ((Delta14C/1000)+1) * exp(-t/lambda)
  if(is.null(er)) {
    return(data.frame(F14C=round(asF, roundby))) 
  } else {
      Fup <- (((Delta14C+er)/1000)+1) * exp(-t/lambda)
      return(round(data.frame(F14C=asF, er=Fup-asF), roundby))
    }
}



#' @name Delta14CtopMC
#' @title Transform Delta14C into pMC
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates Delta14C, the historical level of Delta14C in the year t cal BP, to F14C values. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param Delta14C The Delta14C value to translate
#' @param er Reported error of the Delta14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @param roundby Amount of decimals required for the output. Defaults to \code{roundby=Inf}, no rounding.
#' @return The corresponding F14C value
#' @examples
#'   Delta14CtoF14C(-10, 1, 238)
#' @export
Delta14CtopMC <- function(Delta14C, er=NULL, t, roundby=Inf)
  return(round(Delta14CtoF14C(Delta14C, er, t)*100, roundby))



### functions to be retired below here

#' @name pMC.age
#' @title To be deprecated. Use pMCtoC14 instead.
#' @description Will be deprecated. Use pMCtoC14 instead.
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from pMC values. The reverse function is C14.pMC.
#' @param mn Reported mean of the pMC.
#' @param sdev Reported error of the pMC.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, use \code{1} here.
#' @param decimals Amount of decimals required for the radiocarbon age.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return Radiocarbon ages from pMC values. If pMC values are above 100\%, the resulting radiocarbon ages will be negative.
#' @export
pMC.age <- function(mn, sdev=c(), ratio=100, decimals=0, lambda=8033) {
  message("pMC.age will be deprecated. Use pMCtoC14 instead") 
  pMCtoC14(mn, sdev, decimals, lambda)
}



#' @name age.pMC
#' @title To be deprecated. Use C14topMC instead.
#' @description Calculate pMC values from radiocarbon ages
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#' this function can be used to calculate pMC values from radiocarbon ages. The reverse function of pMC.C14.
#' @param mn Reported mean of the 14C age.
#' @param sdev Reported error of the 14C age.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, a warning is provided; use \code{age.F14C}.
#' @param decimals Amount of decimals required for the pMC value. Defaults to 5.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return pMC values from C14 ages.
#' @export
age.pMC <- function(mn, sdev=c(), ratio=100, decimals=5, lambda=8033) {
  message("age.pMC will be deprecated. Use C14topMC instead") 	
  C14topMC(mn, sdev, decimals, lambda)
}



#' @name F14C.age
#' @title To be deprecated. Calculate C14 ages from F14C values.
#' @description Calculate C14 ages from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C or fraction modern carbon. Since Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param mn Reported mean of the F14C
#' @param sdev Reported error of the F14C. Returns just the mean if left empty.
#' @param decimals Amount of decimals required for the radiocarbon age. Quite sensitive, defaults to 5.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return Radiocarbon ages from F14C values. If F14C values are above 100\%, the resulting radiocarbon ages will be negative.
#' @export
F14C.age <- function(mn, sdev=c(), decimals=5, lambda=8033) {
  message("F14C.age will be deprecated. Use F14CtoC14 instead")
  F14CtoC14(mn, sdev, decimals, lambda)
}



#' @name age.F14C
#' @title To be deprecated. Use C14.F14C instead
#' @description Calculate F14C values from radiocarbon ages
#' @details Post-bomb dates are often reported as F14C or fraction modern carbon. Since Bacon expects radiocarbon ages,
#' this function can be used to calculate F14C values from radiocarbon ages. The reverse function of \link{F14CtoC14}.
#' @param mn Reported mean of the 14C age.
#' @param sdev Reported error of the 14C age. If left empty, will translate mn to F14C.
#' @param decimals Amount of decimals required for the F14C value. Defaults to 5.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return F14C values from C14 ages.
#' @export
age.F14C <- function(mn, sdev=c(), decimals=5, lambda=8033) {
  message("age.F14C will be deprecated. Use C14toF14C instead") 	
  C14toF14C(mn, sdev, decimals, lambda)
}

