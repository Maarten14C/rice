#

# add sample weight functions (per Philippa Ascough's suggestion). Given a %C (perhaps provide estimates for sample types such as peat, bone, ...), a loss during pretreatment, and a required graphite weight, what sample weight will be required?)

# prepare a function to redo deltaR calcs when new Marine curves come out. Using BCADtocalBP(shells$collected), calBPto14C(cc=2) and shells$C14, shells$er. Unclear how the dR errors are obtained.

# fruits-type model that mixes atmospheric and marine calibration curves. Freshwater effects can cause C14 shifts of up to 1k.

# error multipliers, rounding. Could add procedures for different labs, e.g. QUB_bg, etc. This would be useful for reasons of transparency and community standards. Add data from historical UBA standards/backgrounds?



#' @name map.dates
#' @title A map of 180k archaeological C-14 dates
#' @description Produce an interactive, browseable map of the c. 180k (!) archaeological radiocarbon dates provided by the p3k14c R package.
#' @details The p3k14c R package is available on github ("people3k/p3k14c"). See Bocinsky, R. Kyle, Darcy Bird, Lux Miranda, and Jacob Freeman (2022.6). Compendium of R code and data for p3k14c: A synthetic global database of archaeological radiocarbon dates. Accessed 1 July 2025. https://doi.org/10.5281/zenodo.6633635. The p3k14c data has been provided under a CC-0 license.
#' Bird, D., Miranda, L., Vander Linden, M. et al. p3k14c, a synthetic global database of archaeological radiocarbon dates. Sci Data 9, 27 (2022). https://doi.org/10.1038/s41597-022-01118-7
#' You are recommended to have a look at the github instructions for the p3k14c package. The `map.dates` function is provided here for an initial visual analysis of the data; p3k14c provides functions to get much more out of the data.
#' Plotting requires the 'leaflet' R package to be installed. An error will occur if it isn't. If your version of R is below 4.0.0 (released back in April 2020), you're prompted to install a more recent version.
#' @return An interactive, browseable map returning all radiocarbon dates (age and Lab ID) in the p3k14c database.
#' @param S The southern limit of the initial map.
#' @param W The western limit of the initial map.
#' @param N The northern limit of the initial map.
#' @param E The eastern limit of initial map.
#' @param download Whether or not to try to download the p13k14c data. If set to FALSE (default option) and if the p3k14c package is not installed, instructions to do so are provided. If TRUE and p3k14c is not installed, the relevant file will be downloaded from `www.p3k14c.org`. Note: it weighs around 26 MB (c. 4 MB zipped).
#' @param rainbow Whether or not to use a rainbow scale to plot the variable.
#' @param mincol Colour for minimum values. Defaults to 'yellow'.
#' @param maxcol Colour for maximum values. Defaults to 'red'. 
#' @param size Size of the symbols. Defaults to 1.5.
#' @param legend.loc Location of the legend. Defaults to the top right.
#' ## Not run: 
#'   alldates <- map.dates()
# ## End(Not run)
#' @export
map.dates <- function(S=48, W=-15, N=62, E=5, download=FALSE, rainbow=FALSE, mincol="yellow", maxcol="red", size=1.5, legend.loc="topright") {

  # if p3k14c is installed, take the data from there. If not, suggest to install it, or download data from p3k14c.org
  p3k14c <- requireNamespace("p3k14c", quietly=TRUE)
  remotes <- requireNamespace("remotes", quietly=TRUE)
  if(getRversion() < "4.0.0")
    stop("R version 4.0.0 or higher is required for this function. Please upgrade R.")
  lflt <- requireNamespace("leaflet", quietly=TRUE)
  if(!lflt)
    stop("Please install the leaflet package:\ninstall.packages(\"leaflet\")")

  if(p3k14c)
    dates <- p3k14c::p3k14c_data else
      if(download) { # then download from the p3k14c.org site and store in rice's data directory
        dir <- tools::R_user_dir("rice", "data")
        if(!dir.exists(dir))
          dir.create(dir, recursive=TRUE)
        local_file <- file.path(dir, "p3k14c_2022.06.csv")
        url <- "https://www.p3k14c.org/data/p3k14c_2022.06.csv"
        download.file(url, destfile = local_file, mode = "wb")
        dates <- utils::read.csv(local_file)
      } else {
        if(remotes)
          stop("please install p3k14c first: 'remotes::install_github(\"people3k/p3k14c\")', then re-run 'map.dates()'") else
            stop("please install remotes: 'install(\"remotes\")' \nthen install p3k14c: 'remotes::install_github(\"people3k/p3k14c\")' \n then re-run 'map.dates()'")
        }

  dates <- dates[!is.na(dates$Lat),]
  dates <- dates[!is.na(dates$Long),]  
  
  age <- dates$Age
  bins <- cut(age, breaks = 100)
  if(rainbow)
    color_scale <- rainbow(100) else
      color_scale <- grDevices::colorRampPalette(c("yellow", "red"))(100)
  cols <- color_scale[as.numeric(bins)]
  qtiles <- round(quantile(age, probs = c(1, 0.75, 0.5, 0.25, 0)), 1)
  
  hover_labels <- paste0(age, " &plusmn; ", dates$Error, "<br>", dates$Material, "<br>", dates$LabID)

  map <- leaflet::leaflet(data=dates)
  map <- leaflet::fitBounds(map, lng1=W, lat1=S, lng2=E, lat2=N)
  map <- leaflet::addProviderTiles(map, leaflet::providers$Esri.WorldImagery, group = "Esri Satellite")
  map <- leaflet::addCircleMarkers(map, lng=~Long, lat=~Lat, 
    color=cols, radius=size, fillOpacity=0.3, label=lapply(hover_labels, htmltools::HTML))
  map <- leaflet::addLegend(map, position = legend.loc,
    colors=color_scale[c(100, 75, 50, 25, 1)], labels=qtiles, title="C-14 age", opacity = 0.7)
  print(map)
  
  invisible(dates)
}



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
      F <- C14toF14C(age)
  }} else
      F <- C14toF14C(age) # then t is on the C14 scale

  F <- as.numeric(F)
  atoms <- (wght/1e3)*Av/12 # number of C atoms in a mg
  C14 <- as.numeric(round(F * C14.ratio * atoms, 0)) # C14 atoms roundest to nearest number
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



#' @name adjust.fractionation
#' @title Adjust a radiocarbon age for fractionation
#' @description Calculate the radiocarbon age by adjusting a sample's d13C to the reference d13C of -25 permil. It is planned to update this function to more properly reflect calculations in the 14CHRONO lab.
#' @details Radiocarbon ages are corrected for fractionation (which can take place in the field, or during lab pretreatment and measurement), by calculating the radiocarbon age as if the d13C fractionation were at the d13C of the standard (-25 permil). Errors are not taken into account.
#' @return The fractionation-adjusted age.
#' @param y The age of the sample (in C14 by default, but can also be in F or pMC).
#' @param d13C The measured d13C value.
#' @param reference_d13C The reference/standard d13C value (OX2, oxalic acid 2, NIST SRM 4990C made from 1977 French beet molasses), set at -25 permil by default.
#' @param realm Type of radiocarbon age. Can be in `C14` (default), `F14C` or `pMC`.
#' @author Maarten Blaauw
#' @examples
#'   adjust.fractionation(5000, -17)
#' @export
adjust.fractionation <- function(y, d13C, reference_d13C=-25, realm="C14") {
  ratio = (1 + reference_d13C / 1000) / (1 + d13C / 1000)
  realm <- tolower(realm)

  if(grepl("^c", realm))
    return(-8033 * log( ratio^2 * exp(-y / 8033) ))
  if(grepl("^f", realm))
    return(y * ratio^2)
  if(grepl("^p", realm))
    return((y/100) * ratio^2 * 100) # calculate as F

  stop("Unknown realm; use 'C14', 'F', or 'pMC'")
}



#' @name adjust.background
#' @title Adjust a radiocarbon age for background measurements
#' @description Calculate the radiocarbon age by adjusting it for a measured background. It is planned to update this function to more properly reflect calculations in the 14CHRONO lab.
#' @details Radiocarbon ages are measured using a series of standards and backgrounds, and the raw values are then corrected for these background values. Backgrounds are >0 (in F14C) owing to contamination in even the cleanest lab.
#' @return The background-adjusted age.
#' @param y The age of the sample (in C14 by default, but can also be in F or pMC).
#' @param er The error of the date.
#' @param bg The background measurement. Should be in the same realm as that of the sample.
#' @param bg.er The error of the background measurement. Should be in the same realm as that of the sample.
#' @param realm Type of radiocarbon age. Can be in `C14` (default), `F14C` or `pMC`.
#' @author Maarten Blaauw
#' @examples
#'   adjust.background(9000, 50, 45000, 200)
#' @export
adjust.background <- function(y, er, bg, bg.er, realm="C14") {
  realm <- tolower(realm)
  is_c14 <- grepl("^c", realm)
  is_pmc <- grepl("^p", realm)
  is_f14c <- grepl("^f", realm)

  if((is_c14 && y > bg) || ((is_f14c || is_pmc) && y < bg))
    stop("sample's age is older than background age!")

  if(is_c14) {
    tmp <- C14toF14C(y, er)
    y <- tmp[,1]; er <- tmp[,2]
    tmp <- C14toF14C(bg, bg.er)
    bg <- tmp[,1]; bg.er <- tmp[,2]
  } else
    if(is_pmc) {
      tmp <- pMCtoF14C(y, er)
      y <- tmp[,1]; er <- tmp[,2]
      tmp <- pMCtoF14C(bg, bg.er)
      bg <- tmp[,1]; bg.er <- tmp[,2]
    } else
    if(!is_f14c)
      stop("Unknown realm; use 'C14', 'F', or 'pMC'")

  X <- y - bg
  D <- 1 - bg
  sigma_X <- sqrt(er^2 + bg.er^2)

  F_corr <- X / D
  er_corr <- sqrt( (1/D)^2 * sigma_X^2 + (X/(D^2))^2 * bg.er^2 )

  if(is_f14c)
    return(data.frame(F14C=F_corr, F14C_er=er_corr)) else
    if(is_pmc)
      return(data.frame(pMC=100*F_corr, pMC_er=100*er_corr)) else
      if(is_c14) {
        as.C <- data.frame(F14CtoC14(F_corr, er_corr))
        return(data.frame(C14=as.C[,1], C14.er=as.C[,2]))
      }
}


