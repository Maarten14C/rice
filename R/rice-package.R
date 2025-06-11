#' @name rice
#' @title Radiocarbon Equations
#' @keywords internal
"_PACKAGE"

#' @docType package
#' @aliases rice-package
"_PACKAGE"
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk>
#' @description Provides equations for the handling of radiocarbon dates, for example to calculate different radiocarbon realms (C14 age, F14C, pMC, D14C), for their calibration, and estimating the effects of contamination. This package accompanies the data package rintcal.
#' @importFrom utils read.table write.table packageName data installed.packages browseURL
#' @importFrom stats approx dnorm median weighted.mean runif pchisq density quantile dgamma rgamma rnorm rbeta sd
#' @importFrom grDevices rgb extendrange grey rainbow colorRampPalette 
#' @importFrom graphics axis par legend lines points polygon segments text mtext abline image rect curve arrows strwidth layout
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot geom_sf coord_sf geom_point aes scale_color_gradient scale_color_gradientn labs theme element_line element_rect
#' @importFrom maps map
#' @importFrom rintcal ccurve intcal.data new.ccdir list.ccurves glue.ccurves mix.ccurves
#' @importFrom ggplot2 geom_sf ggplot coord_sf theme element_line element_rect 
#' @name rice

## usethis namespace: start
## usethis namespace: end
NULL

# trying to deal with reported NOTE:
#   "Found the following files/directories: ‘rnaturalearthhires’"
.onLoad <- function(libname, pkgname) {
  if (is.null(Sys.getenv("R_RNATURAL_EARTH_CACHE", unset = NA))) {
    Sys.setenv("R_RNATURAL_EARTH_CACHE" = tempdir())
  }
}

#
# # function from rintcal, which for rice's fromto function requires as.D:
# glue.ccurves <- function(prebomb="IntCal20", postbomb="NH1", thisprebombcurve=c(), thispostbombcurve=c(), as.F=FALSE, as.pMC=FALSE, as.D=FALSE, cc.dir=c(), decimals=8) {
#   if(length(thispostbombcurve) == 0)
#     postbomb <- ccurve(postbomb, TRUE, cc.dir=cc.dir, as.F=as.F, as.pMC=as.pMC, as.D=as.D, decimals=decimals) else
#       postbomb <- thispostbombcurve
#   if(length(thisprebombcurve) == 0)
#     prebomb <- ccurve(prebomb, FALSE, cc.dir=cc.dir, as.F=as.F, as.pMC=as.pMC, as.D=as.D, decimals=decimals) else
#       prebomb <- thisprebombcurve
#
#   glued <- rbind(postbomb, prebomb)
#   glued <- glued[order(glued[,1]),]
#   repeated <- which(diff(glued[,1]) == 0)
#   if(length(repeated) > 0)
#     invisible(glued[-repeated,]) else # remove any repeated years
#       invisible(glued[order(glued[,1]),])
# }
