#' @name rice
#' @title Radiocarbon Calibration Equations
#' @keywords internal
"_PACKAGE"

#' @docType package
#' @aliases rice-package
"_PACKAGE"
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk>
#' @description Provides equations for the calibration of radiocarbon dates, as well as to calculate different radiocarbon realms (C14 age, F14C, pMC, D14C) and estimating the effects of contamination. This package accompanies the data package rintcal.

#' @importFrom utils read.table write.table packageName
#' @importFrom stats approx dnorm median weighted.mean
#' @importFrom grDevices rgb extendrange grey rainbow
#' @importFrom graphics axis par legend lines points polygon segments text mtext abline image
#' @importFrom data.table fread fwrite
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom rintcal ccurve intcal.data intcal.read.data intcal.write.data intcal.data.frames list.ccurves mix.ccurves glue.ccurves new.ccdir
#' @name rice

## usethis namespace: start
## usethis namespace: end
NULL
