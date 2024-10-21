
# from/to	calBP		BCAD		C14			F14C		pMC			D14C	
# calBP					calBPtoBCAD	calBPtoC14	calBPtoF14C	calBPtopMC	calBPtoD14C
# BCAD		BCADtocalP				BCADtoC14	BCADtoF14C	BCADtopMC	BCADtoD14C
# C14		NA			NA						C14toF14C	C14topMC	C14toD14C
# F14C		NA			NA			F14CtoC14				F14CtopMC	F14CtoD14C
# pMC		NA			NA			pMCtoC14	pMCtoF14C				pMCtoD14C
# D14C		NA			NA			D14CtoC14	D14CtoF14C	D14CtopMC



#' @name calBPtoBCAD
#' @title calculate BC/AD ages from cal BP ages
#' @details Turn cal BP ages into BC/AD (or BCE/CE). Negative ages indicate BC, positive ages AD. Since the Gregorian and Julian calendars do not include 0 BCAD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from 1 (AD) to -1 (i.e., 1 BC). Other calendars, such as the astronomical one, do include zero. The often-used BCE/CE ages are equivalent to BC/AD.
#' @param x The calBP age(s) to be translated into BC/AD ages. 
#' @param zero Whether or not zero BC/AD should be included. Defaults to 0. 
#' @return The BC/AD age(s). BC ages are negative, AD ages are positive.
#' @examples
#'  calBPtoBCAD(2024)
#'  calBPtoBCAD(1945:1955, zero=TRUE)
#'  calBPtoBCAD(1945:1955, zero=FALSE)
#' @export 
calBPtoBCAD <- function(x, zero=TRUE) {
  if(zero) 
    return(1950-x) else {
      X <- 1950-x 
      neg <- which(x >= 1950)  
      if(length(neg) > 0)
        X[neg] <- 1949 - x[neg]
      return(X) # x = 0 is not an error
    }
}



#' @name calBPtoC14
#' @title Find the 14C age and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding 14C age and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param yr The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#' calBPtoC14(100)
#' @export
calBPtoC14 <- function(yr, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL) {
  if(is.null(thiscurve))
    cc <- rintcal::ccurve(cc=cc, postbomb=postbomb, cc.dir=cc.dir) else
      cc <- thiscurve
  mu <- approx(cc[,1], cc[,2], yr, rule=rule)$y
  er <- approx(cc[,1], cc[,3], yr, rule=rule)$y
  return(cbind(mu, er, deparse.level=0))
}



#' @name calBPtoF14C
#' @title Find the F14C and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param yr The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   calBPtoF14C(100)
#' @export
calBPtoF14C <- function(yr, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL) {
  y <- calBPtoC14(yr, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve)
  return(C14toF14C(y[,1], y[,2]))
}



#' @name calBPtopMC
#' @title Find the pMC and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param yr The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   calBPtopMC(100)
#' @export
calBPtopMC <- function(yr, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL) {
  y <- calBPtoC14(yr, cc, postbomb, rule, cc.dir, thiscurve)
  return(C14topMC(y[,1], y[,2]))
}



#' @name calBPtoD14C
#' @title Find the pMC and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param yr The cal BP year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   calBPtoD14C(100)
#' @export
calBPtoD14C <- function(yr, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL) {
  F <- calBPtoF14C(yr, cc, postbomb, rule, cc.dir, thiscurve)
  Dmn <- F14CtoD14C(F[,1], t=yr)
  Dup <- F14CtoD14C(F[,1]+F[,2], t=yr)
  return(cbind(Dmn, Dup-Dmn, deparse.level=0))
}



#' @name BCADtocalBP
#' @title calculate cal BP ages from BC/AD ages
#' @details Turn BC/AD (or BCE/CE) ages into cal BP ages. Negative ages indicate BC, positive ages AD. Since the Gregorian and Julian calendars do not include 0 BC/AD (i.e., 31 December of 1 BC is followed by 1 January of AD 1), zero can be omitted. The years then go from 1 (AD) to -1 (i.e., 1 BC). Other calendars, such as the astronomical one, do include zero. The often-used BCE/CE ages are equivalent to BC/AD.
#' @param x The BCAD age(s) to be translated into cal BP age(s). BC ages are negative, AD ages are positive.
#' @param zero Whether or not zero BC/AD should be included. Defaults to 0. 
#' @return The cal BP age(s).
#' @examples
#'  BCADtocalBP(2024)
#'  BCADtocalBP(-1, zero=TRUE)
#'  BCADtocalBP(-1, zero=FALSE)
#' @export 
BCADtocalBP <- function(x, zero=TRUE) 
  if(zero)
    return(1950-x) else {
      X <- 1950-x 
      neg <- which(x < 0)  
      if(length(neg) > 0)
        X[neg] <- 1949-x[neg]
        return(X) # x = 0 is not an error
    }



#' @name BCADtoC14
#' @title Find the 14C age and error belonging to a BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding 14C age and error are returned. BC ages are negative. In this implementation, the year 0 BC/AD does exist.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than AD 1950, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered BC/AD age
#' @param yr The BC/AD year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param zero Whether or not to include 0 in BC/AD years. Defaults to TRUE.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   BCADtoC14(100)
#' @export
BCADtoC14 <- function(yr, cc=1, postbomb=FALSE, zero=TRUE, rule=1, cc.dir=NULL, thiscurve=NULL) 
  return(calBPtoC14(BCADtocalBP(yr, zero=zero), 
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve))



#' @name BCADtoF14C
#' @title Find the F14C and error belonging to a BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C and error are returned. BC ages are negative. In this implementation, the year 0 BC/AD does exist.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than AD 1950, a postbomb curve will have to be provided.
#' @return The calibration-curve F14C belonging to the entered BC/AD age
#' @param yr The BC/AD year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param zero Whether or not to include 0 in BC/AD years. Defaults to TRUE.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   BCADtoF14C(100)
#' @export
BCADtoF14C <- function(yr, cc=1, postbomb=FALSE, zero=TRUE, rule=1, cc.dir=NULL, thiscurve=NULL) 
  return(calBPtoF14C(BCADtocalBP(yr, zero=zero), 
    cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve))



#' @name BCADtopMC
#' @title Find the pMC and error belonging to a BC/AD age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding pMC and error are returned. BC ages are negative. In this implementation, the year 0 BC/AD does exist.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For ages younger than AD 1950, a postbomb curve will have to be provided.
#' @return The calibration-curve F14C belonging to the entered BC/AD age
#' @param yr The BC/AD year.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param zero Whether or not to include 0 in BC/AD years. Defaults to TRUE.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   BCADtopMC(100)
#' @export
BCADtopMC <- function(yr, cc=1, postbomb=FALSE, zero=TRUE, rule=1, cc.dir=NULL, thiscurve=NULL) 
    return(calBPtopMC(BCADtocalBP(yr, zero=zero), 
      cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve))



#' @name BCADtoD14C
#' @title Find the pMC and error belonging to a cal BP age.
#' @description Given a calendar age, the calibration curve (default cc=1) is interpolated and the corresponding F14C value and error are returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided.
#' @return The calibration-curve 14C year belonging to the entered cal BP age
#' @param yr The cal BP year.
#' @param zero Whether or not to include 0 in BC/AD years. Defaults to TRUE.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). 
#' @author Maarten Blaauw
#' @examples
#'   BCADtoD14C(1900)
#' @export
BCADtoD14C <- function(yr, zero=TRUE, cc=1, postbomb=FALSE, rule=1, cc.dir=NULL, thiscurve=NULL) {
  calBP <- BCADtocalBP(yr, zero)
  Fres <- calBPtoF14C(calBP, cc=cc, postbomb=postbomb, rule=rule, cc.dir=cc.dir, thiscurve=thiscurve)
  Dmn <- F14CtoD14C(Fres[,1], t=yr)
  Dup <- F14CtoD14C(Fres[,1]+Fres[,2], t=yr)
  return(cbind(Dmn, Dup-Dmn, deparse.level=0))
}



#' @name C14toF14C
#' @title Calculate F14C values from C14 ages
#' @description Calculate F14C values from radiocarbon ages
#' @details Post-bomb dates are often reported as F14C or fraction modern carbon. Since software such as Bacon expects radiocarbon ages,
#' this function can be used to calculate F14C values from radiocarbon ages. The reverse function of \link{F14C.age}.
#' @param mn Reported mean of the 14C age.
#' @param sdev Reported error of the 14C age. If left empty, will translate mn to F14C.
#' @param decimals Amount of decimals required for the F14C value. Defaults to 5.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return F14C values from C14 ages.
#' @examples
#'   C14toF14C(-2000, 20)
#' @export
C14toF14C <- function(mn, sdev=NULL, decimals=5, lambda=8033) {
  y <- exp(-mn / lambda)
  if(is.null(sdev))
    return(signif(y, decimals)) else {
      sdev <- y - exp(-(mn + sdev) / lambda)
      return(signif(cbind(y, sdev, deparse.level=0), decimals))
    }
}



#' @name C14topMC
#' @title Calculate pMC values from C14 ages
#' @description Calculate pMC values from radiocarbon ages
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#' this function can be used to calculate pMC values from radiocarbon ages. The reverse function of \link{pMCtoC14}.
#' @param mn Reported mean of the C14 age.
#' @param sdev Reported error of the C14 age.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, a warning is provided; use \code{C14.F14C}.
#' @param decimals Amount of decimals required for the pMC value. Defaults to 5.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return pMC values from C14 ages.
#' @examples
#'   C14topMC(-2000, 20)
#'   C14topMC(-2000, 20, 1)
#' @export
C14topMC <- function(mn, sdev=NULL, ratio=100, decimals=5, lambda=8033) {
  if(ratio !=100)
    warning("C14topMC expects a ratio of 100. For ratio=1, use C14toF14C")
  y <- exp(-mn / lambda)
  if(is.null(sdev))
    return(signif(ratio*y, decimals)) else {
      sdev <- y - exp(-(mn + sdev) / lambda)
      return(signif(ratio*cbind(y, sdev, deparse.level=0), decimals))
  }
}



#' @name C14toD14C
#' @title Transform C14 age(s) into D14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates C14 ages into Delta14C, the historical level of Delta14C in the year t cal BP. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param C14 The C14 age to translate
#' @param sdev Reported error of the C14 age. Returns just the mean if left empty.
#' @param t the cal BP age
#' @return The corresponding D14C value
#' @examples
#'   C14toD14C(0.985, 20, 222)
#' @export
C14toD14C <- function(C14, sdev=NULL, t) {
  asF <- cbind(C14toF14C(C14, sdev))
  Dmn <- 1000 * ((asF[,1] / exp(-t/8267)) - 1)
  if(is.null(sdev))
    return(Dmn) else {
	  Fup <- asF[,1] + asF[,2]	
	  Dup <- 1000 * ((Fup / exp(-t/8267)) - 1)
	  return(cbind(Dmn, Dup-Dmn, deparse.level=0))
    }     
}



#' @name F14CtoC14
#' @title Calculate C14 ages from F14C values.
#' @description Calculate C14 ages from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param mn Reported mean of the F14C
#' @param sdev Reported error of the F14C. Returns just the mean if left empty.
#' @param decimals Amount of decimals required for the radiocarbon age. Quite sensitive, defaults to 5.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return The radiocarbon ages from the F14C values. If F14C values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   F14CtoC14(1.10, 0.5) # a postbomb date, so with a negative C14 age
#'   F14CtoC14(.80, 0.5) # prebomb dates can also be calculated
#' @export
F14CtoC14 <- function(mn, sdev=NULL, decimals=5, lambda=8033) {
  y <- -lambda * log(mn)
  if(is.null(sdev))
    signif(y, decimals) else {
    sdev <- y - -lambda * log((mn+sdev))
    signif(cbind(y, sdev, deparse.level=0), decimals)
  }
}



#' @name F14CtopMC
#' @title Calculate pMC ages from F14C values.
#' @description Calculate pMC values from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param mn Reported mean of the F14C
#' @param sdev Reported error of the F14C. Returns just the mean if left empty.
#' @return The pMC values from the F14C values. Basically the original values multiplied by 100.
#' @examples
#'   F14CtopMC(1.10, 0.5)
#' @export
F14CtopMC <- function(mn, sdev=NULL)
  return(100*cbind(mn, sdev, deparse.level=0))



#' @name F14CtoD14C
#' @title Transform F14C into D14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates F14C values into Delta14C, the historical level of Delta14C in the year t cal BP. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param F14C The F14C value to translate
#' @param sdev Reported error of the F14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @return The corresponding D14C value
#' @examples
#'   F14CtoD14C(0.89, .001, 900)
#' @export
F14CtoD14C <- function(F14C, sdev=NULL, t) {
  Dmn <- 1000 * ((F14C / exp(-t/8267)) - 1)
  if(is.null(sdev))
    return(Dmn) else {
      Dup <- 1000 * (((F14C+sdev) / exp(-t/8267)) - 1)
	  return(cbind(Dmn, Dup-Dmn, deparse.level=0))
    }
}



#' @name pMCtoC14
#' @title Calculate C14 ages from pMC values.
#' @description Calculate C14 ages from pMC values of radiocarbon dates.
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from pMC values. The reverse function is C14.pMC.
#' @param mn Reported mean of the pMC.
#' @param sdev Reported error of the pMC.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, use \code{1} here.
#' @param decimals Amount of decimals required for the radiocarbon age.
#' @param lambda The mean-life of radiocarbon (based on Libby half-life of 5568 years)
#' @return Radiocarbon ages from pMC values. If pMC values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   pMCtoC14(110, 0.5) # a postbomb date, so with a negative 14C age
#'   pMCtoC14(80, 0.5) # prebomb dates can also be calculated
#'   pMCtoC14(.8, 0.005, ratio=1) # throws a warning, use F14C.age instead
#' @export
pMCtoC14 <- function(mn, sdev=NULL, ratio=100, decimals=0, lambda=8033) { 
  if(ratio !=100)
    warning("pMCtoC14 expects a ratio of 100. For ratio=1, use F14CtoC14")
  y <- -lambda * log(mn/ratio)
  if(is.null(sdev))
    signif(y, decimals) else {
    sdev <- y - -lambda * log((mn+sdev)/ratio)
    round(cbind(y, sdev, deparse.level=0), decimals)
  }
}



#' @name pMCtoF14C
#' @title Calculate pMC ages from F14C values.
#' @description Calculate pMC values from F14C values of radiocarbon dates.
#' @details Post-bomb dates are often reported as F14C (between 0 at c. 55 kcal BP and 1 at c. AD 1950). Since software such as Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from F14C values. The reverse function is \link{age.F14C}.
#' @param mn Reported mean of the F14C
#' @param sdev Reported error of the pMC value. Returns just the mean if left empty.
#' @return The F14C values from the pMC values. Basically the original values divided by 100.
#' @examples
#'   pMCtoF14C(110, 5)
#' @export
pMCtoF14C <- function(mn, sdev=NULL)
  return(cbind(mn, sdev, deparse.level=0)/100)



#' @name pMCtoD14C
#' @title Transform F14C into D14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates F14C values into Delta14C, the historical level of Delta14C in the year t cal BP. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param pMC The pMC value to translate
#' @param sdev Reported error of the pMC value. Returns just the mean if left empty.
#' @param t the cal BP age
#' @return The corresponding D14C value
#' @examples
#'   pMCtoD14C(0.985, .1, 222)
#' @export
pMCtoD14C <- function(pMC, sdev=NULL, t) {
  asF <- pMCtoF14C(pMC, sdev)	
  F14CtoD14C(asF[,1], asF[,2], t)
}


#' @name D14CtoC14
#' @title Transform D14C into C14 age
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates Delta14C, the historical level of Delta14C in the year t cal BP, to C14 ages. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param D14C The Delta14C value to translate
#' @param sdev Reported error of the D14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @return The corresponding C14 age
#' @examples
#'   D14CtoC14(-10, 1, 238)
#' @export
D14CtoC14 <- function(D14C, sdev=NULL, t) {
  return( F14CtoC14( D14CtoF14C(D14C=D14C, sdev=sdev, t=t) ))
}



#' @name D14CtoF14C
#' @title Transform D14C into F14C
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates Delta14C, the historical level of Delta14C in the year t cal BP, to F14C values. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param D14C The Delta14C value to translate
#' @param sdev Reported error of the D14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @return The corresponding F14C value
#' @examples
#'   D14CtoF14C(-10, 1, 238)
#' @export
D14CtoF14C <- function(D14C, sdev=NULL, t) {
  asF <- ((D14C/1000)+1) * exp(-t/8267)
  if(is.null(sdev))
    return(asF) else {
      Fup <- (((D14C+sdev)/1000)+1) * exp(-t/8267)
	  return(cbind(asF, Fup-asF, deparse.level=0))
    }
}



#' @name D14CtopMC
#' @title Transform D14C into pMC
#' @details As explained by Heaton et al. 2020 (Radiocarbon), 14C measurements are commonly expressed in
#' three domains: Delta14C, F14C and the radiocarbon age. This function translates Delta14C, the historical level of Delta14C in the year t cal BP, to F14C values. Note that per convention, this function uses the Cambridge half-life, not the Libby half-life.
#' @param D14C The Delta14C value to translate
#' @param sdev Reported error of the D14C. Returns just the mean if left empty.
#' @param t the cal BP age
#' @return The corresponding F14C value
#' @examples
#'   D14CtoF14C(-10, 1, 238)
#' @export
D14CtopMC <- function(D14C, sdev=NULL, t) 
  return(D14CtoF14C(D14C, sdev, t)/100)



### functions to be retired below here

#' @name pMC.age
#' @title To be deprecated. Use pMCtoC14 instead.
#' @description Will soon be deprecated. Use pMCtoC14 instead.
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
  pMCtoC14(mn=mn, sdev=sdev, ratio=ratio, decimals=decimals, lambda=lambda)
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
  C14topMC(mn=mn, sdev=sdev, ratio=ratio, decimals=decimals, lambda=lambda)
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
  F14CtoC14(mn=mn, sdev=sdev, decimals=decimals, lambda=lambda)
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
  C14toF14C(mn=mn, sdev=sdev, decimals=decimals, lambda=lambda)
}

