#' @name rice
#' @title Radiocarbon Equations
#' @keywords internal
"_PACKAGE"

#' @docType package
#' @aliases rice-package
"_PACKAGE"
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk>
#' @description Provides equations for the handling of radiocarbon dates, for example to calculate different radiocarbon realms (C14 age, F14C, pMC, D14C), for their calibration, and estimating the effects of contamination. This package accompanies the data package rintcal.

#' @importFrom utils read.table write.table packageName data installed.packages
#' @importFrom stats approx dnorm median weighted.mean runif
#' @importFrom grDevices rgb extendrange grey rainbow
#' @importFrom graphics axis par legend lines points polygon segments text mtext abline image rect
#' @importFrom rlang sym
#' @importFrom rnaturalearth ne_countries
#' @importFrom ggplot2 ggplot geom_sf coord_sf geom_point aes scale_color_gradient scale_color_gradientn labs theme element_line element_rect
#' @name rice

## usethis namespace: start
## usethis namespace: end
NULL

require(rintcal)