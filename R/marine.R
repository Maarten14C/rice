ocean.map <- function(S, W, N, E, scale = c(), ocean.col = "aliceblue", land.col = rgb(0, 0.5, 0, 0.6)) {
  # Check if required libraries are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required")
  }
  
  if(!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
    message("For high-resolution maps, install rnaturalearthhires from GitHub using: devtools::install_github('ropensci/rnaturalearthhires').")
  }  
  
  ihave <- installed.packages()
  if(("rnaturalearth" %in% ihave) && ("sf" %in% ihave)) {
    # Handle different scales and package availability
    if (length(scale) == 0 || scale == "medium") {
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    } else if (scale == "small") {
      world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
    } else if (scale == "large") {
      if (!"rnaturalearthhires" %in% installed.packages()) {
        message("Using a medium-scale map. For higher resolution, install rnaturalearthhires (devtools::install_github(\"ropensci/rnaturalearthhires\")).")
        world <- rnaturalearthdata::countries50
      } else {
        world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
      }
    }

    p <- ggplot(data = world) +
           geom_sf(fill = land.col) +
           coord_sf(xlim = c(W, E), ylim = c(S, N), expand = TRUE) +
           theme(
             panel.grid.major = element_line(color = rgb(0, 0, 0, 0.5), linetype = 2, size = 0.1),
             panel.background = element_rect(fill = ocean.col),
             legend.background = element_rect(fill = "transparent"),
             legend.key = element_rect(fill = "transparent")
           )
  } else {
      message("rnaturalearth is not installed - plotting a basic map.")
	  message("Please issue the command install.packages('rnaturalearth'), then try finding/mapping shells again")
	  p <- c()
    }
  return(p)	
}



#' @name find.shells
#' @title Find nearby shell-derived dR values
#' @description Find the shells closest to a chosen coordinate, and plot the dR values and feeding ecology. Uses the marine database downloaded (30 Aug 2024) from calib.org/marine. See Reimer PJ, Reimer RW, 2001. A marine reservoir correction database and on-line interface. Radiocarbon 43:461-3.
#' @return A dataset with the n nearest dR values, and a plot of their coordinates.
#' @param longitude Longitude of the point. Can only deal with one point at a time.
#' @param latitude Latitude of the point. Can only deal with one point at a time.
#' @param nearest The number of shell values to be returned. Defaults to 50.
#' @param colour The variable to be plotted as colour. Expects a continuous variable. Defaults to 'dR'.
#' @param rainbow Whether or not to use a rainbow scale to plot the variable.
#' @param size Size of the symbols. Defaults to 2.
#' @param scale Resolution of the map. Can be "small", "medium" or "large". If the latter, a high-resolution dataset will have to be downloaded using the R package 'rnaturalearthhires'. Since this package is on github but not on CRAN, you will have to download it yourself (using the command devtools::install_github("ropensci/rnaturalearthhires")). Defaults to 'medium' if 'rnaturalearthhires' is not installed, and to 'high' if it is installed.
#' @param mincol Colour for minimum values.
#' @param maxcol Colour for maximum values.
#' @param symbol The variable to be plotted as symbol. Expects a categoric variable. Defaults to 'feeding'.
#' @param symbol.legend Whether or not to plot the legend for the symbols.
#' @param ocean.col Colour for the oceans. Defaults to \code{ocean.col="aliceblue"}.
#' @param land.col Colour for the land. Defaults to semi-transparent darkgreen: \code{land.col=rgb(0, 0.5, 0, 0.6)}.
#' @examples
#'   UK <- find.shells(0, 55, scale="medium")
#'   mean(UK$dR)
#'   Caribbean <- find.shells(-70, 20, 30, scale="medium")
#' @export
find.shells <- function(longitude, latitude, nearest=50, colour='dR', rainbow=FALSE, size=2, scale=c(), mincol="yellow", maxcol="red", symbol='feeding', symbol.legend=TRUE, ocean.col="aliceblue", land.col=rgb(0, 0.5, 0., 0.6)) {
  lon <- lat <- NULL # to get rid of subsequent ggplot2-related warnings
  if(length(c(longitude,latitude)) != 2)
    stop("we need 1 entry for longitude, 1 for latitude")

  shells <- get("shells", envir = .GlobalEnv)

  # from https://stackoverflow.com/questions/27928/calculate-distance-between-two-latitude-longitude-points-haversine-formula
  hav.dist <- function(long1, lat1, long2, lat2) {
    R <- 6371
    p <- pi/180
    diff.long <- (long2 - long1) * p
    diff.lat <- (lat2 - lat1) * p
    a <- sin(diff.lat/2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin(diff.long/2)^2
    b <- 2 * asin(pmin(1, sqrt(a)))
    d = R * b
    return(d)
  }

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

  xl <- abs(W-E)/20
  yl <- abs(N-S)/20
  xmid <- longitude
  ymid <- latitude

  p <- ocean.map(S, W, N, E, scale, ocean.col, land.col)
  
  if(length(p) == 0) { # then rnaturalearth is not installed, so we plot a basic map instead
	padding <- 1 
  	maps::map(xlim=c(W-padding, E+padding), ylim=c(S-padding, N+padding), fill = TRUE, col = land.col, bg = ocean.col)
	color_scale <- grDevices::colorRampPalette(c("yellow", "red"))(100)
	cols <- color_scale[as.numeric(cut(nearshells[,5], breaks = 100))]
	points(nearshells[,1], nearshells[,2], col=cols, pch=20)
  } else {

    if(symbol.legend)
      p <- p +
        geom_point(data=nearshells, aes(x=lon, y=lat, color=!!sym(colour), shape=!!sym(symbol)), size=size, alpha=.8) else
          p <- p +
            geom_point(data=nearshells, aes(x=lon, y=lat, color=!!sym(colour)), size=size, alpha=.8)

    if(rainbow)
      p <- p + scale_color_gradientn(colors = rainbow(7)) + labs(shape="feeding") else
        p <- p + scale_color_gradient(low=mincol, high=maxcol) + labs(shape="feeding")

    print(p)
  }

  return(nearshells)
}



#' @name map.shells
#' @title Plot regional shell-derived dR values
#' @description Find the shells that fit within a rectangular region (bounded by N, E, S and W), and plot the dR values and feeding ecology. Uses the marine database downloaded (30 Aug 2024) from calib.org/marine. See Reimer PJ, Reimer RW, 2001. A marine reservoir correction database and on-line interface. Radiocarbon 43:461-3. Expects the coordinates for the map to be provided (starting south, then clockwise as with R axes).
#' @return A plot and the relevant dR values.
#' @param S The southern limit of the rectangular region.
#' @param W The western limit of the rectangular region.
#' @param N The northern limit of the rectangular region.
#' @param E The eastern limit of the rectangular region.
#' @param colour The variable to be plotted as colour. Expects a continuous variable. Defaults to 'dR'.
#' @param rainbow Whether or not to use a rainbow scale to plot the variable.
#' @param size Size of the symbols. Defaults to 2.
#' @param scale Resolution of the map. Can be "small", "medium" or "large". If the latter, a high-resolution dataset will have to be downloaded using the R package 'rnaturalearthhires'. Since this package is on github but not on CRAN, you will have to download it yourself (using the command devtools::install_github("ropensci/rnaturalearthhires")). Defaults to 'medium' if 'rnaturalearthhires' is not installed, and to 'high' if it is installed.
#' @param mincol Colour for minimum values.
#' @param maxcol Colour for maximum values.
#' @param symbol The variable to be plotted as symbol. Expects a categoric variable. Defaults to 'feeding'. 
#' @param symbol.legend Whether or not to plot the legend for the symbols.
#' @param ocean.col Colour for the oceans. Defaults to \code{ocean.col="aliceblue"}.
#' @param land.col Colour for the land. Defaults to semi-transparent darkgreen: \code{land.col=rgb(0, 0.5, 0, 0.6)}.
#' @examples
#'  N_UK <- map.shells(53, -11, 60, 2, scale="medium")
#'  mean(N_UK$dR)
#' @export
map.shells <- function(S=48,W=-15, N=62, E=5, colour='dR', rainbow=FALSE, size=2, scale=c(), mincol="yellow", maxcol="red", symbol='feeding', symbol.legend=TRUE, ocean.col="aliceblue", land.col=rgb(0, 0.5, 0., 0.6)) {
  lon <- lat <- NULL # to get rid of subsequent ggplot2-related warnings
  shells <- get("shells", envir = .GlobalEnv)
  shells[[symbol]] <- as.factor(shells[[symbol]])
  sel <- shells[shells$lon>=W & shells$lon<=E & shells$lat>=S & shells$lat <= N,]

  p <- ocean.map(S, W, N, E, scale, ocean.col, land.col)

  if(length(p) == 0) { # then rnaturalearth is not installed, so we plot a basic map instead
	padding <- 1 
  	maps::map(xlim=c(W-padding, E+padding), ylim=c(S-padding, N+padding), fill = TRUE, col = land.col, bg = ocean.col)
	color_scale <- grDevices::colorRampPalette(c("yellow", "red"))(100)
	cols <- color_scale[as.numeric(cut(sel[,5], breaks = 100))]  # Assign colors based on the value
	points(sel[,1], sel[,2], col=cols, pch=20)
  } else {

  if(symbol.legend) 
    p <- p +
      geom_point(data=sel, aes(x=lon, y=lat, color=!!sym(colour), shape=!!sym(symbol)), size=size, alpha=.8) else
        p <- p +
          geom_point(data=sel, aes(x=lon, y=lat, color=!!sym(colour)), size=size, alpha=.8)

  if(rainbow)
    p <- p + scale_color_gradientn(colors = rainbow(7)) + labs(shape="feeding") else
      p <- p + scale_color_gradient(low=mincol, high=maxcol) + labs(shape="feeding")

    print(p)
  } 
  
  return(sel)
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
#'   N_UK <- map.shells(53, -11, 60, 2, scale="medium")
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
#' @examples
#'  N_UK <- map.shells(53, -11, 60, 2, scale="medium")
#'  shells.mean(N_UK)
#'  nearby <- find.shells(0,56,20) # somewhere in Scotland
#'  shells.mean(nearby, distance=TRUE) # distance matters
#' @export
shells.mean <- function(dat, feeding=c(), draw=TRUE, distance=FALSE, pch=20, col.mn=1, lty.mn=2, col.sd=rgb(0,0,0,.1)) {
  if(length(feeding)>0)
    dat <- dat[which(dat$feeding %in% feeding)]
  wmn_sd <- weighted_means(dat$dR, dat$dSTD)
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
