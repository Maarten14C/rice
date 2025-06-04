# add b2k to the realms equations

# titles of contaminate plots were plotted outside of the range - now using par(xpd=TRUE)

# can ocean circulation model output be added to ocean maps? A la https://svs.gsfc.nasa.gov/vis/a000000/a003800/a003821/flat_ocean07we_4096rgb.0001-05_print.jpg?
# or https://earth.nullschool.net/#current/ocean/surface/currents/overlay=significant_wave_height/orthographic=7.26,70.02,3106/loc=-19.763,50.473

# check how is.pMC and is.F work in calibrate(). Make interpolation to e.g. years more intelligent (default c() then 1 if prebomb, .1 if postbomb

# rintcal has as.F through ccurve(as.F=TRUE)) (but not as.pMC)

# add data from historical UBA standards/backgrounds?

# terr-marine contribution calculation

# AMS background and fractionation corrections



#' @name howmanyC14
#' @title Amount of C14 particles in a sample
#' @description Find the amount of remaining C14 atoms in a sample, given its weight and age.
#' @details The number of carbon atoms in the sample is estimated. Given the known C14/C ratio at F=1, and given the sample's age, we can estimate the number of remaining C14 atoms.
#' @return The estimated number of C14 atoms.
#' @param age The age of the sample (in cal BP per default, or in C14 BP is use.cc=FALSE).
#' @param wght The weight of the sample (in mg). Defaults to 1 mg.
#' @param use.cc Whether or not to use the calibration curve. If set to \code{use.cc=FALSE}, then we assume that the age is the radiocarbon age (this enables ages beyond the reach of the calibration curves to be used).
#' @param Av Avogadro's number, used to calculate the number of carbon atoms in the sample.
#' @param C14.ratio The 14C/C ratio at F=1 (AD 1950).
#' @param format The format of the printed numbers. Defaults to either scientific (for large numbers) or as fixed-point, depending on the size of the number.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param talk Whether or not to provide feedback (defaults to TRUE).
#' @param decimals Number of decimals to be returned for F and atom counts.
#' @author Maarten Blaauw
#' @examples
#'   howmanyC14(0) # recent sample
#'   howmanyC14(55e3) # at dating limit
#'   howmanyC14(145e3) # way beyond the dating limit, 1 C14 atom per mg remains
#' @export
howmanyC14 <- function(age, wght=1, use.cc=TRUE, Av=6.02214076e23, C14.ratio=1.176e-12, format="g", cc=1, postbomb=FALSE, cc.dir=NULL, thiscurve=NULL, talk=TRUE, decimals=3) {

  if(use.cc) {
    F <- calBPtoF14C(age, cc=cc, postbomb=postbomb, cc.dir=cc.dir, thiscurve=thiscurve)[,1]
    if(is.na(F)) {
      message("Cannot use calibration curve for this age, assuming C14 age")
      F <- C14toF14C(age, decimals=10)
  }} else
      F <- C14toF14C(age, decimals=10) # then t is on the C14 scale

  atoms <- (wght/1e3)*Av/12 # number of C atoms in a mg
  C14 <- round(F * C14.ratio * atoms, 0) # C14 atoms roundest to nearest number
  perminute <- round(C14/wght/30,0)
  persecond <- round(perminute/60,0)

  atoms <- formatC(atoms, format=format, digits=decimals)
  C14.talk <- formatC(C14, format=format, digits=decimals)
  
  decays <- round(C14 * log(2) / (5730 * 365.25), decimals)
  decays <- formatC(decays, format=format, digits=decimals)

  if(talk) {
    message(wght, " mg carbon contains ", atoms, " C atoms")
    message("C14 atoms remaining at ", age, " cal BP (F=", round(F, decimals), "): ", C14.talk)
	message(decays, " C-14 atoms in the sample will decay each day")
    message("For a 1 mg AMS target, assuming a 100% efficiency, ", perminute, " particles would be counted per minute, or ", persecond, " per second")
  }

  invisible(C14)
}

