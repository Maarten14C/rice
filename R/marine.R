
ocean.map <- function(S, W, N, E, shells=c(), browse=FALSE, mapsize="large", padding=0.1, ocean.col="aliceblue", land.col=rgb(0, 0.5, 0, 0.6), rainbow=FALSE, symbol='feeding', symbol.legend=TRUE, legend.loc=c(.95, .02), legend.size=c(.05, .20), mincol="yellow", maxcol="red", colour='dR', warn=TRUE) {

  # Check if useful libraries are installed
  hassf <- requireNamespace("sf", quietly=TRUE)
  rne <- requireNamespace("rnaturalearth", quietly=TRUE)
  rnedata <- requireNamespace("rnaturalearthdata", quietly=TRUE)
  remotes <- requireNamespace("remotes", quietly=TRUE)
  lflt <- requireNamespace("leaflet", quietly=TRUE)
  coper <- requireNamespace("CopernicusMarine", quietly=TRUE)
  hiresmaps <- requireNamespace("rnaturalearthhires", quietly = TRUE)

  if(!hiresmaps && warn && mapsize == "large") {
    message("High-resolution maps not available; using standard resolution.")
	mapsize <- "medium"
  }

  if(warn) {
    if(browse) {
      if(getRversion() < "4.0.0")
        stop("For browseable maps, R >= 4.0.0 is required. Please update your R installation.")		
      if(!lflt)
        stop("Please install the leaflet package:\ninstall.packages(\"leaflet\")")
      if(!coper)
        stop("Please install the CopernicusMarine package:\ninstall.packages(\"CopernicusMarine\")")
    } else {
        if(!hassf)
          message("Using (ugly) basic map. For better maps, please install the packages sf and rnaturalearth: \ninstall.packages(c(\"sf\", \"rnaturalearth\"))") else {
          if(mapsize=="large") {
            if(!rne)
              message("Please install the rnaturalearth package:",
                "\ninstall.packages(\"rnaturalearth\")")
            if(!rnedata)
              message("Please install the rnaturalearthdata package:",
                "\ninstall.packages(\"rnaturalearthdata\")")

            # rnaturalearthhires is nice but has to be installed from github
            if(!hiresmaps)
              if(remotes)
                message("For detailed maps, install rnaturalearthhires from GitHub:\n",
                  "`remotes::install_github('ropensci/rnaturalearthhires')`") else
                    message("Install first remotes and then rnaturalearthhires:\n",
                      "install.packages(\"remotes\")\n",
                      "`remotes::install_github(\"ropensci/rnaturalearthhires\")`")
          }
        }
      }
  }

  if(rainbow)
    color_scale <- rainbow(100) else
      color_scale <- grDevices::colorRampPalette(c("yellow", "red"))(100)

  if(browse) {
    dR <- shells[,5]
    bins <- cut(dR, breaks = 100)
    cols <- color_scale[as.numeric(bins)]
    qtiles <- round(quantile(dR, probs = c(1, 0.75, 0.5, 0.25, 0)), 1)

    hover_labels <- paste0("&Delta;R: ", dR, " &plusmn; ", shells[,6],
      "<br>No.: ", shells$no, "<br>Feeding: ", shells$feeding)
    map <- leaflet::leaflet()
    map <- leaflet::addProviderTiles(map, leaflet::providers$Esri.WorldImagery, group = "Esri Satellite")
    map <- leaflet::addLayersControl(map,
      baseGroups = c("Esri Satellite", "Positron", "Sea Water Potential Temperature"),
      options = leaflet::layersControlOptions(collapsed = FALSE))
    map <- CopernicusMarine::addCmsWMTSTiles(map,
      product = "GLOBAL_ANALYSISFORECAST_PHY_001_024",
      layer = "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m",
      variable = "thetao",
      tilematrixset = "EPSG:3857",
      options = leaflet::WMSTileOptions(format = "image/png", transparent = TRUE),
      group = "Sea Water Potential Temperature"
    )
    map <- leaflet::addLayersControl(map,
      baseGroups = c("ESRI Satellite", "Sea Water Potential Temperature"),
      options = leaflet::layersControlOptions(collapsed = FALSE)
    )
    map <- leaflet::addCircleMarkers(map, data=shells, lng=~lon, lat=~lat,
      color=cols, radius=5, fillOpacity=0.8,
      label=lapply(hover_labels, htmltools::HTML))
    map <- leaflet::addLegend(map, position = "bottomright",
      colors=color_scale[c(100, 75, 50, 25, 1)],
      labels=qtiles, title="dR", opacity = 0.7)
    print(map)
    return() # don't continue on to following code
  }

  if(!hassf) { # then we'll plot basic maps (ugly)
    width <- abs(E-W); height <- abs(N-S)
    par(mar=rep(1,4))
    maps::map(xlim=c(W-(padding*width), E+(padding*width)), ylim=c(S-(padding*height), N+(padding*height)),
      fill = TRUE, col = land.col, bg = ocean.col)
    cols <- color_scale[as.numeric(cut(shells[,5], breaks = 100))]
    points(shells[,1], shells[,2], col=cols, pch=20)

    value_range <- range(shells[, 5], na.rm = TRUE)

    coors <- par("usr")
    xmin <- coors[1] + legend.loc[1] * (coors[2] - coors[1])
    xmax <- coors[1] + (legend.loc[1]+legend.size[1]) * (coors[2] - coors[1])
    ymin <- coors[3] + legend.loc[2] * (coors[4] - coors[3])
    ymax <- coors[3] + (legend.loc[2]+legend.size[2]) * (coors[4] - coors[3])

    yticks <- seq(ymin, ymax, length.out=4)
    vals <- round(seq(min(shells[,5]), max(shells[,5]), length=length(yticks)), 0)
    image(x=c(xmin, xmax), y=seq(ymin, ymax, length.out=100), z=matrix(1:100, nrow = 1),
      col=color_scale, axes=FALSE, add=TRUE)
    text((xmin+xmax)/2, ymax, expression(Delta*R), adj=c(0.5,0))
    text(xmax+.2, yticks, vals, cex=.5, adj=c(0, .5))

    return() # plotted basic map, end of function
  }

  if (!hassf || !rne) {
    world <- sf::st_as_sf(maps::map("world", fill=TRUE, plot=FALSE))
  } else {
    if(mapsize=="large") # this is not working as expected, find.shells(0, 55, mapsize="large")
      if(!hiresmaps)
        mapsize <- "medium"
    world <- rnaturalearth::ne_countries(scale=mapsize, returnclass="sf")
  } 

  p <- ggplot(data = world) +
    geom_sf(fill = land.col) +
    coord_sf(xlim = c(W, E), ylim = c(S, N), expand = TRUE) +
      theme(
        panel.grid.major = element_line(color = rgb(0, 0, 0, 0.5), linetype = 2, linewidth = 0.1),
        panel.background = element_rect(fill = ocean.col),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")
      )

  lon_col <- sym("lon")
  lat_col <- sym("lat")

  if(symbol.legend)
    p <- p + geom_point(data=shells, aes(x=!!lon_col, y=!!lat_col, color=!!sym(colour), shape=!!sym(symbol)), size=2, alpha=.8) else
      p <- p + geom_point(data=shells, aes(x=!!lon_col, y=!!lat_col, color=!!sym(colour)), size=2, alpha=.8)

  if(rainbow) {
    p <- p + scale_color_gradientn(colors = rainbow(7)) + labs(color = expression(Delta*R), shape = "feeding")
  } else
      p <- p + scale_color_gradient(low = mincol, high = maxcol) + labs(color = expression(Delta*R), shape = "feeding")

  print(p)
}



# from https://stackoverflow.com/questions/27928/calculate-distance-between-two-latitude-longitude-points-haversine-formula
hav.dist <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth's diameter in km
  p <- pi/180
  diff.long <- (long2 - long1) * p
  diff.lat <- (lat2 - lat1) * p
  a <- sin(diff.lat/2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin(diff.long/2)^2
  b <- 2 * asin(pmin(1, sqrt(a)))
  return(R * b)
}



#' @name find.shells
#' @title Find nearby shell-derived dR values
#' @description Find the shells closest to a chosen coordinate, and plot the dR values and feeding ecology. Uses the marine database downloaded (30 Aug 2024) from calib.org/marine. See Reimer PJ, Reimer RW, 2001. A marine reservoir correction database and on-line interface. Radiocarbon 43:461-3.
#' @details
#' This function uses the `rnaturalearth` package for country maps. If the high-resolution maps are desired,
#' the `rnaturalearthhires` package must be installed from GitHub.
#' @return A dataset with the n nearest dR values, and a plot of their coordinates.
#' @param longitude Longitude of the point. Can only deal with one point at a time.
#' @param latitude Latitude of the point. Can only deal with one point at a time.
#' @param nearest The number of shell values to be returned. Defaults to 50.
#' @param browse Type of map to provide. \code{browse=FALSE} (default) plots a static map in R's device (doesn't require Internet access), while \code{browse=TRUE} opens a browsable, interactive map in your Internet browser.
#' @param colour The variable to be plotted as colour. Expects a continuous variable. Defaults to 'dR'.
#' @param rainbow Whether or not to use a rainbow scale to plot the variable.
#' @param size Size of the symbols. Defaults to 2.
#' @param mapsize Resolution of the map. Can be "small" or "large". If the latter, a high-resolution dataset will have to be downloaded using the R package 'rnaturalearthhires'. Since this package is on github but not on CRAN, you will have to download it yourself (using the command remotes::install_github("ropensci/rnaturalearthhires")). Defaults to "small" if 'rnaturalearthhires' is not installed, and to "large" if it is installed.
#' @param mincol Colour for minimum values.
#' @param maxcol Colour for maximum values.
#' @param feeding Optionally, the output of only specific types of feeding ecology (e.g., deposit, suspension, browser) can be selected. Defaults to returning all feeding ecologies.
#' @param symbol The variable to be plotted as symbol. Expects a categoric variable. Defaults to 'feeding'.
#' @param symbol.legend Whether or not to plot the legend for the symbols.
#' @param legend.loc Location of the legend, if using a basic plot. Defaults to the bottom right corner based on par("usr"), \code{legend.loc=c(0.95, 0.02)}
#' @param legend.size Size of the legend, if using a basic plot. Defaults to \code{legend.size=c(0.05, 0.2)}
#' @param ocean.col Colour for the oceans. Defaults to \code{ocean.col="aliceblue"}.
#' @param land.col Colour for the land. Defaults to semi-transparent darkgreen: \code{land.col=rgb(0, 0.5, 0, 0.6)}.
#' @param padding Area around the map if using a basic plot. Avoids strange line features. Defaults to \code{padding=1}. 
#' @param warn Whether or not to warn if some recommended R packages are not available.
#' @param currents If set to TRUE (the default), the user will be asked if they want to browse a map of ocean currents. If the user responds 'y', an Internet browser window will be opened pointing to a zoomed-in map of ocean currents (at 50 m depth). The ocean currents are from 'earth.nullschool.net' and are based on an ocean circulation model which is updated daily. Owing to limitations of the website, the shell locations cannot currently be added to the page itself.
#' @examples
#'   UK <- find.shells(0, 55, mapsize="small")
#'   mean(UK$dR)
#'   Caribbean <- find.shells(-70, 20, 30, mapsize="small")
#' @export
find.shells <- function(longitude, latitude, nearest=50, browse=FALSE, colour="dR", rainbow=FALSE, size=2, mapsize="large", mincol="yellow", maxcol="red", feeding=c(), symbol="feeding", symbol.legend=TRUE, legend.loc=c(0.95, 0.02), legend.size=c(0.05, 0.2), ocean.col="aliceblue", land.col=rgb(0, 0.5, 0., 0.6), padding=1, warn=TRUE, currents=TRUE) {
  lon <- lat <- NULL # to get rid of subsequent ggplot2-related warnings
  if(length(c(longitude,latitude)) != 2)
    stop("we need 1 entry for longitude, 1 for latitude")

  if(!is.logical(currents) || length(currents) != 1)
    stop("'currents' must be a single TRUE or FALSE")
  
  shells <- get("shells", envir = .GlobalEnv)
  if(length(feeding) > 0)
    if(tolower(feeding) %in% c("suspension", "unknown", "carnivore", "browser", "deposit", "scavenger", "commensalism", "omnivore"))
      shells <- shells[shells$feeding == feeding,] else
        warning("cannot find this feeding ecology") 

  shell_coors <- data.frame(lon=shells$lon, lat=shells$lat)
  distances <- apply(shell_coors, 1, function(row) {
    hav.dist(longitude, latitude, row["lon"], row["lat"])})

  o <- order(distances)
  o <- o[1:nearest]
  distances <- distances[o]
  nearshells <- cbind(shells[o,], distances)
  
  S <- min(nearshells$lat)
  W <- min(nearshells$lon)
  N <- max(nearshells$lat)
  E <- max(nearshells$lon)

  ocean.map(S, W, N, E, shells=nearshells,
    mapsize=mapsize, browse=browse, ocean.col=ocean.col, land.col=land.col,
    rainbow=rainbow, symbol=symbol, symbol.legend=symbol.legend, 
	legend.loc=legend.loc, legend.size=legend.size, 
	mincol=mincol, maxcol=maxcol, colour=colour, 
	warn=warn, padding=padding)
  
  if(interactive())
    if(isTRUE(currents)) {
      if(E < W) E <- E + 360
      width <- hav.dist(W, (S+N)/2, E, (S+N)/2)
      height <- hav.dist((W+E)/2, S, (W+E)/2, N)
      w <- 10000 * (1000 / width)
      h <- 10000 * ((1000 / 2) / height) 
      sz <- max(300, min(5000, ceiling(min(w, h))))
      answer <- tolower(readline("Do you want to browse a map with ocean currents? [y/N]: "))
      if(answer == "y") {
        latitude  <- round((S + N) / 2, 3)
        longitude <- round((W + E) / 2, 3)
        browseURL(paste0("https://earth.nullschool.net/#current/ocean/surface/currents/orthographic=",
          round(longitude, 3), ",", round(latitude, 3), ",", sz))
      }
    }  

  invisible(nearshells)
}



#' @name map.shells
#' @title Plot regional shell-derived dR values
#' @description Find the shells that fit within a rectangular region (bounded by N, E, S and W), and plot the dR values and feeding ecology. Uses the marine database downloaded (30 Aug 2024) from calib.org/marine. See Reimer PJ, Reimer RW, 2001. A marine reservoir correction database and on-line interface. Radiocarbon 43:461-3. Expects the coordinates for the map to be provided (starting south, then clockwise as with R axes).
#' @details
#' This function uses the `rnaturalearth` package for country maps. If the high-resolution maps are desired,
#' the `rnaturalearthhires` package must be installed from GitHub.
#' @return A plot and the relevant dR values.
#' @param S The southern limit of the rectangular region.
#' @param W The western limit of the rectangular region.
#' @param N The northern limit of the rectangular region.
#' @param E The eastern limit of the rectangular region.
#' @param browse Type of map to provide. \code{browse=FALSE} (default) plots a static map in R's device (doesn't require Internet access), while \code{browse=TRUE} opens a browsable, interactive map in your Internet browser.
#' @param colour The variable to be plotted as colour. Expects a continuous variable. Defaults to 'dR'.
#' @param rainbow Whether or not to use a rainbow scale to plot the variable.
#' @param size Size of the symbols. Defaults to 2.
#' @param mapsize Resolution of the map. Can be "small" or "large". If the latter, a high-resolution dataset will have to be downloaded using the R package 'rnaturalearthhires'. Since this package is on github but not on CRAN, you will have to download it yourself (using the command remotes::install_github("ropensci/rnaturalearthhires")). Defaults to "small" if 'rnaturalearthhires' is not installed, and to "large" if it is installed.
#' @param mincol Colour for minimum values.
#' @param maxcol Colour for maximum values.
#' @param feeding Optionally, the output of only specific types of feeding ecology (e.g., deposit, suspension, browser) can be selected. Defaults to returning all feeding ecologies.
#' @param symbol The variable to be plotted as symbol. Expects a categoric variable. Defaults to 'feeding'. 
#' @param symbol.legend Whether or not to plot the legend for the symbols.
#' @param ocean.col Colour for the oceans. Defaults to \code{ocean.col="aliceblue"}.
#' @param land.col Colour for the land. Defaults to semi-transparent darkgreen: \code{land.col=rgb(0, 0.5, 0, 0.6)}.
#' @param legend.loc Location of the legend, if using a basic plot. Defaults to the bottom right corner based on par("usr"), \code{legend.loc=c(0.95, 0.02)}
#' @param legend.size Size of the legend, if using a basic plot. Defaults to \code{legend.size=c(0.05, 0.2)}
#' @param padding Area around the map if using a basic plot. Avoids strange line features. Defaults to \code{padding=0.1}. 
#' @param warn Whether or not to warn if some recommended R packages are not available.
#' @param currents If set to TRUE (the default), the user will be asked if they want to browse a map of ocean currents. If the user responds 'y', an Internet browser window will be opened pointing to a zoomed-in map of ocean currents (at 50 m depth). The ocean currents are from 'earth.nullschool.net' and are based on an ocean circulation model which is updated daily. Owing to limitations of the website, the shell locations cannot currently be added to the page itself.
#' @examples
#'  N_UK <- map.shells(53, -11, 60, 2, mapsize="small")
#'  mean(N_UK$dR)
#' @export
map.shells <- function(S=48, W=-15, N=62, E=5, browse=FALSE, colour="dR", rainbow=FALSE, size=2, mapsize="large", mincol="yellow", maxcol="red", feeding=c(), symbol="feeding", symbol.legend=TRUE, ocean.col="aliceblue", land.col=rgb(0, 0.5, 0., 0.6), legend.loc=c(.95, .02), legend.size=c(.05, .2), padding=0.1, warn=TRUE, currents=TRUE) {
  lon <- lat <- NULL # to get rid of subsequent ggplot2-related warnings
  shells <- get("shells", envir = asNamespace("rice")) # envir was .GlobalEnv
  shells[[symbol]] <- as.factor(shells[[symbol]])
  if(length(feeding) > 0)
    if(tolower(feeding) %in% c("suspension", "unknown", "carnivore", "browser", "deposit", "scavenger", "commensalism", "omnivore"))
      shells <- shells[shells$feeding == feeding,] else
        warning("cannot find this feeding ecology")

  sel <- shells[shells$lon>=W & shells$lon<=E & shells$lat>=S & shells$lat <= N,]

  if(!is.logical(currents) || length(currents) != 1)
    stop("'currents' must be a single TRUE or FALSE") 

  ocean.map(S, W, N, E, shells=sel,
    mapsize=mapsize, browse=browse, ocean.col=ocean.col, 
	  land.col=land.col, rainbow=rainbow, symbol=symbol, 
	  symbol.legend=symbol.legend, legend.loc=legend.loc, 
	  legend.size=legend.size, mincol=mincol, maxcol=maxcol, 
	  colour=colour, warn=warn, padding=padding)
  
  if(interactive())
    if(isTRUE(currents)) {
      if(E < W) E <- E + 360
      width <- hav.dist(W, (S+N)/2, E, (S+N)/2)
      height <- hav.dist((W+E)/2, S, (W+E)/2, N)
      w <- 10000 * (1000 / width)
      h <- 10000 * ((1000 / 2) / height) 
      sz <- max(300, min(5000, ceiling(min(w, h))))
      answer <- tolower(readline("Do you want to browse a map with ocean currents? [y/N]: "))
      if(answer == "y") {
        latitude  <- round((S + N) / 2, 3)
        longitude <- round((W + E) / 2, 3)
        browseURL(paste0("https://earth.nullschool.net/#current/ocean/surface/currents/orthographic=",
          round(longitude, 3), ",", round(latitude, 3), ",", sz))
        }
      } 
  
  invisible(sel)
}  



#' @name weighted_means
#' @title Calculate the weighted mean of C14 ages
#' @description Calculating the weighted mean of multiple C14 ages, using their means and lab errors.  
#' @return The weighted mean and error (the latter is the maximum of the weighted error and the square root of the variance).
#' @param y The C14 ages.
#' @param er The lab errors of the C14 ages.
#' @param round Rounding to be applied (defaults to 1 decimal).
#' @param talk Report details of the found values.
#' @examples
#'   N_UK <- map.shells(53, -11, 60, 2, mapsize="small")
#'   weighted_means(N_UK$dR, N_UK$dSTD)
#' @export
weighted_means <- function(y, er, round=1, talk=TRUE) {
  if(any(er == 0))
    stop("errors need to be >0")
  if(sum(is.na(er)) > 0)
    stop("we cannot have NAs for errors") 
  if(length(y) - length(er) != 0)
    stop("we need as many ages as errors")
  wghts <- 1 / (er^2)
  w_mean <- sum(y * wghts) / sum(wghts)
  w_er <- round(sqrt(1 / sum(wghts)), round)
  
  # also calculate the weighted variance of dR
  n_eff <- (sum(wghts))^2 / sum(wghts^2)
  w_var <- sum(wghts * (y - w_mean)^2) / sum(wghts) * (n_eff / (n_eff - 1))
  sdev <- round(sqrt(w_var), round)
  w_mean <- round(w_mean, round)
  
  if(talk)
    message(paste0("wmean ", w_mean, ", error is max of weighted uncertainty (", w_er, ") & sdev (", sdev, "): ", max(w_er, sdev)))
  return(c(w_mean, max(w_er, sdev)))
}



#' @name shells.mean
#' @title Plot and summarize the dR values
#' @description After selecting a relevant range of shell values, plot them and calculate the weighted mean and variance.  
#' @return A plot of the dR values, as well as the weighted mean (vertical line) and (weighted) error (rectangle).
#' @param dat The data, as returned from the function 'plot.shells'.
#' @param feeding Whether or not to select a specific feeding behaviour. Defaults to empty (no selection of feeding behaviour).
#' @param draw Whether or not to draw the values. 
#' @param distance Plot the dR values according to their distance (if you've used find.shells; assumes that 'dat' has a final column with the distances).
#' @param pch Symbol to be plotted. Defaults to a closed circle (\code{pch=20}).
#' @param col.mn Colour for the weighted mean. Defaults to black, \code{col.mn=1}.
#' @param lty.mn Line type for the weighted mean. Defaults to dashed, \code{lty.mn=2}.
#' @param col.sd Colour of the rectangle of the error. Defaults to transparent grey, \code{col.sd=rgb(0,0,0,.1)}.
#' @param talk Report details of the found values.
#' @examples
#'  N_UK <- map.shells(53, -11, 60, 2, mapsize="small")
#'  shells.mean(N_UK)
#'  nearby <- find.shells(0,56,20) # somewhere in Scotland
#'  shells.mean(nearby, distance=TRUE) # distance matters
#' @export
shells.mean <- function(dat, feeding=c(), draw=TRUE, distance=FALSE, pch=20, col.mn=1, lty.mn=2, col.sd=rgb(0,0,0,.1), talk=TRUE) {
  if(length(feeding)>0)
    dat <- dat[which(dat$feeding %in% feeding)]
  wmn_sd <- weighted_means(dat$dR, dat$dSTD, talk=talk)
  if(draw) {
    food <- dat$feeding
    allfood <- unique(food)
    dye <- c()
    for(i in 1:length(allfood))
      dye[which(food==allfood[i])] <- i
    rng <- range(dat$dR-dat$dSTD, dat$dR+dat$dSTD)
    if(distance)
      yvals <- dat[,ncol(dat)] else
        yvals <- 1:nrow(dat)
    if(distance)
      ylab <- "distance (km)" else ylab=""
    plot(dat$dR, yvals, pch=pch, col=dye, xlim=rng, xlab="delta R", ylab=ylab)
    segments(dat$dR-dat$dSTD, yvals, dat$dR+dat$dSTD, yvals, col=dye)
    abline(v=wmn_sd[1], lty=lty.mn, col=col.mn)
    rect(wmn_sd[1]-wmn_sd[2], min(yvals)-10, wmn_sd[1]+wmn_sd[2], max(yvals)+10, col=col.sd, border=NA)
    legend("topleft", legend=allfood, text.col=1:length(allfood), bty="n", cex=.5)
  }

  return(wmn_sd)
} 
