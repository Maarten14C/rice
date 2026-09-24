# how to get the latest version of the shells data is described in data-raw/extractfrommarinedatabase.R

#' shells Data
#'
#' A dataset containing the deltaR values and accompanying data from the marine database
#'
#' @docType data
#' @format A data frame with 1968 rows and 15 variables.
#' \describe{
#'   \item{lon}{Longitude of the datapoint}
#'   \item{lat}{Latitude of the datapoint}
#'   \item{no}{Map or ID number of the datapoint}
#'   \item{taxonN}{Taxon number of the datapoint}
#'   \item{dR}{calculated deltaR of the datapoint}
#'   \item{dSTD}{uncertainty of the deltaR of the datapoint}
#'   \item{collected}{Collection year for the datapoint}
#'   \item{res}{Reservoir effect of the datapoint}
#'   \item{res.error}{Uncertainty of the reservoir effect of the datapoint}
#'   \item{C14}{Radiocarbon age of the datapoint}
#'   \item{er}{Error of the radiocarbon age of the datapoint}
#'   \item{lab}{Lab code of the datapoint}
#'   \item{ref}{Reference for the datapoint}
#'   \item{taxon}{Taxon of the datapoint}
#'   \item{feeding}{Feeding ecology of the datapoint (if known)}
#' }
#' @source Data downloaded from calib.org/marine
#' @examples
#' data(shells)
#' head(shells)
"shells"



# elements data
# # downloaded from https://www-nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all
# list_elements <- read.csv("data-raw/IAEA_APIoutput.csv")
# list_elements <- list_elements[-(1:3),] # first three rows are not necessary
# elementnames <- read.csv("data-raw/elements.csv")
# abundance <- list_elements$abundance
# stable <- which(list_elements$half_life == "STABLE")
# protons <- list_elements$z
# neutrons <- list_elements$n
# radiocarbon <- which(protons==6 & neutrons==8)
# abundance[radiocarbon] <- 1.2e-12 # correcting the values for C14!
# sym <- list_elements$symbol
# el.name <- c()
# for(i in 1:length(sym))
#   el.name[i] <- elementnames[which(elementnames[,2] == sym[i]),3]
# halfl <- list_elements$half_life_sec
# halfl[which(list_elements$half_life == "STABLE")] <- 0
# halfl[radiocarbon] <- 1.80825e+11 # correct halflife of C14 (in seconds)
# decay <- list_elements$decay_1
# decay[decay==""] <- NA
# decay[stable] <- "stable"

# elements <- list(
#   protons   = as.numeric(protons),
#   neutrons  = as.numeric(neutrons),
#   name      = as.character(el.name),
#   symbol    = as.character(sym),
#   decay     = as.character(decay),
#   halflife  = as.numeric(halfl),
#   abundance = as.numeric(abundance)
# )
# save(elements, file="data/elements.rda", compress="bzip2")

#' elements Data
#'
#' A dataset containing the amounts of protons and neutrons, names, symbols, decay type and halflife for all naturally occurring isotopes
#'
#' @docType data
#' @format A list with 7 variables and 282 entries for each variable.
#' \describe{
#'   \item{protons}{amount of protons in the nucleus}
#'   \item{neutrons}{amount of neutrons in the nucleus}
#'   \item{name}{name of the element}
#'   \item{symbol}{symbol of the element}
#'   \item{decay}{type of decay product (if not stable)}
#'   \item{halflife}{decay half-life}
#'   \item{abundance}{abundance of the isotope, relative to the element}
#' }
#' @source Data downloaded from https://www-nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all
#' @examples
#' data(elements)
#' head(elements)
"elements"


