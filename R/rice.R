# to do: 

# check how only relevant functions from the rintcal package can be loaded. This also because rice will be loaded by clam, rbacon, rplum, coffee

# add data from historical UBA standards/backgrounds?

# soil F14C simulation from Vegard

# lifetime function that models F_tot based on C input and removal for each year, given known t_start, t_end, F_t

# terr-marine contribution calculation

# AMS background and fractionation corrections

# add 'realm' to calibrate function? Perhaps better not, since then the date would also have to be on that realm

lifetime <- function(birth, death, add.mean, add.shape, recycle.mean, recycle.shape, recyle.max=0.05, cc=1, cc.dir=c(), thiscurve=c(), by=1, BCAD=FALSE, visualise=TRUE) {

  if(BCAD) {
    if(death <= birth)
      stop("birth year must precede death year")
    birth <- BCADtocalBP(birth)
  } else
      if(death >= birth)
        stop("birth year must precede death year")
  yr <- seq(0, death-birth, by=by)
  F.lifetime <- calBPtoF14C(yr, cc, cc.dir=cc.dir, thiscurve=thiscurve)

  cbind(yr, F.lifetime) # how deal with the cc errors?


}
