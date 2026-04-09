

#' @name map.dates
#' @title A map of 180k archaeological C-14 dates
#' @description Produce an interactive, browseable map of the c. 180k (!) archaeological radiocarbon dates provided by the p3k14c R package.
#' @details The p3k14c R package is not on CRAN but is available on github ("people3k/p3k14c"). See Bocinsky, R. Kyle, Darcy Bird, Lux Miranda, and Jacob Freeman (2022.6). Compendium of R code and data for p3k14c: A synthetic global database of archaeological radiocarbon dates.  https://doi.org/10.5281/zenodo.6633635. The p3k14c data has been provided under a CC-0 license.
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
  remotes <- requireNamespace("remotes", quietly=TRUE)
  if(getRversion() < "4.0.0")
    stop("R version 4.0.0 or higher is required for this function. Please upgrade R.")
  lflt <- requireNamespace("leaflet", quietly=TRUE)
  if(!lflt)
    stop("Please install the leaflet package:\ninstall.packages(\"leaflet\")")

  # fedora-gcc on github does not like requireNamespace rnaturalearthhires
  p3k14c <- nzchar(system.file(package="p3k14c"))

  if(p3k14c)
    dates <- get("p3k14c_data", envir=asNamespace("p3k14c")) else
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
