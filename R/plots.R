#' @name draw.ccurve
#' @title Draw a calibration curve.
#' @description Draw one or two of the calibration curves, or add a calibration curve to an existing plot.
#' @return A plot of the calibration curve
#' @param cal1 First calendar year for the plot. Defaults to 0 cal BP.
#' @param cal2 Last calendar year for the plot. Defaults to 55,000 cal BP.
#' @param cc1 Name of the calibration curve. Can be "IntCal20", "Marine20", "SHCal20", or for the previous curves "IntCal13", "Marine13" or "SHCal13". Can also be "nh1", "nh2", "nh3", "sh1-2", "sh3", "nh1_monthly", "nh1_monthly", "nh2_monthly", "nh3_monthly", "sh1-2_monthly", "sh3_monthly", "Kure", "LevinKromer" or "Santos" for postbomb curves.
#' @param cc2 Optional second calibration curve to plot. Can be "IntCal20", "Marine20", "SHCal20", or for the previous curves "IntCal13", "Marine13" or "SHCal13". Defaults to nothing, NA.
#' @param cc1.postbomb Use \code{postbomb=TRUE} to get a postbomb calibration curve for cc1 (default \code{cc1.postbomb=FALSE}).
#' @param cc2.postbomb Use \code{postbomb=TRUE} to get a postbomb calibration curve for cc2 (default \code{cc2.postbomb=FALSE}).
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param realm Which 'realm' of radiocarbon to use. Defaults to \code{realm="C14"} but can also be set to \code{realm="F14C"}, \code{realm="pMC"} or \code{realm="D14C"}. Can be shorted to, respectively, "C", "F", "P" or "D" (or their lower-case equivalents).
#' @param realm2 Which 'realm' to use for the second calibration curve (if used). Defaults to \code{realm="C14"} but can also be set to \code{realm="F14C"}, \code{realm="pMC"} or \code{realm="D14C"}. Can be shorted to, respectively, "C", "F", "P" or "D" (or their lower-case equivalents).
#' @param cal.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}), or to \code{age.lab="kcal BP"} etc. if ka=TRUE.
#' @param cal.rev Reverse the calendar axis. 
#' @param c14.lab Label for the C-14 axis. Defaults to 14C BP (or 14C kBP if ka=TRUE).
#' @param c14.lim Axis limits for the C-14 axis. Calculated automatically by default. 
#' @param c14.rev Reverse the C-14 axis.
#' @param ka Use kcal BP (and C14 kBP).
#' @param add.yaxis Whether or not to plot the second calibration. Defaults to \code{add.yaxis=FALSE}.
#' @param cc1.col Colour of the calibration curve (outline).
#' @param cc1.fill Colour of the calibration curve (fill).
#' @param cc2.col Colour of the calibration curve (outline), if activated (default cc2=NA).
#' @param cc2.fill Colour of the calibration curve (fill), if activated (default cc2=NA).
#' @param add Whether or not to add the curve(s) to an existing plot. Defaults to FALSE, which draws a new plot
#' @param bty Draw a box around a box of a certain shape. Defaults to \code{bty="l"}.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param legend Location of the legend (only activated if more than one curve is plotted). Plotted in the topleft corner by default. Use \code{legend=c()} to leave empty
#' @param ... Any additional optional plotting parameters. 
#' @examples 
#' draw.ccurve()
#' draw.ccurve(1000, 3000, cc2="Marine20")
#' draw.ccurve(1800, 2020, BCAD=TRUE, cc2="nh1", cc2.postbomb=TRUE)
#' draw.ccurve(1800, 2010, BCAD=TRUE, cc2="nh1", add.yaxis=TRUE)
#' @export
draw.ccurve <- function(cal1=c(), cal2=c(), cc1="IntCal20", cc2=NA, cc1.postbomb=FALSE, cc2.postbomb=FALSE, BCAD=FALSE, realm="C14", realm2=c(), cal.lab=NA, cal.rev=FALSE, c14.lab=NA, c14.lim=NA, c14.rev=FALSE, ka=FALSE, add.yaxis=FALSE, cc1.col=rgb(0,0,1,.5), cc1.fill=rgb(0,0,1,.2), cc2.col=rgb(0,.5,0,.5), cc2.fill=rgb(0,.5,0,.2), add=FALSE, bty="l", cc.dir=NULL, legend="topleft", ...) {

  # read and narrow down the calibration curve(s)
  if(cc1 %in% c(2, "Marine20")) # then no postbomb curve available
    cc.1 <- rintcal::ccurve(2, postbomb=FALSE, cc.dir) else
      if(cc1.postbomb)
        cc.1 <- rintcal::glue.ccurves(cc1, cc1.postbomb, cc.dir) else
          cc.1 <- rintcal::ccurve(cc1, cc1.postbomb, cc.dir)
  cc.cal <- 1 # which column to use for calendar ages
  if(BCAD) {
    cc.1[,4] <- 1950 - cc.1[,1] # add a column...
    cc.cal <- 4 # ... and use it
  }
  if(length(cal1) == 0)
    cal1 <- cc.1[1,cc.cal]
  if(length(cal2) == 0)
    cal2 <- cc.1[nrow(cc.1),cc.cal]
  
  mindat <- cc.1[,cc.cal] >= min(cal1, cal2)
  maxdat <- cc.1[,cc.cal] <= max(cal1, cal2)
  cc.1 <- cc.1[which(mindat * maxdat == 1),]

  if(ka) {
    cc.1[,1] <- cc.1[,1]/1e3
    if(grepl("c", tolower(realm)))  # ka doesn't make sense for F14C, pMC, or D14C
      cc.1[,2:3] <- cc.1[,2:3]/1e3
  }

  if(grepl("^f", tolower(realm))) {
    F <- C14toF14C(cc.1[,2], cc.1[,3])
    cc.1[,2:3] <- F
  }
  if(grepl("^p", tolower(realm))) {
    p <- C14topMC(cc.1[,2], cc.1[,3])
    cc.1[,2:3] <- p
  }
  if(grepl("^d", tolower(realm))) {
    asD <- C14toD14C(cc.1[,2], cc.1[,3], cc.1[,1])
    cc.1[,2:3] <- cbind(asD)
  }
 
  cc1.pol <- cbind(c(cc.1[,cc.cal], rev(cc.1[,cc.cal])), c(cc.1[,2]-cc.1[,3], rev(cc.1[,2]+cc.1[,3])))
  
  if(!is.na(cc2)) {
    if(length(realm2) == 0)
      realm2 <- realm
    if(cc2.postbomb)
      cc.2 <- rintcal::glue.ccurves(cc2, cc2.postbomb, cc.dir) else
        cc.2 <- rintcal::ccurve(cc2, cc2.postbomb, cc.dir)
    if(BCAD) 
      cc.2[,4] <- 1950 - cc.2[,1] 
    mindat <- cc.2[,cc.cal] >= .9*min(cal1, cal2)
    maxdat <- cc.2[,cc.cal] <= 1.1*max(cal1, cal2)

    if(ka) {
      cc.2[,1] <- cc.2[,1]/1e3
      if(grepl("c", tolower(realm2)))  # ka doesn't make sense for F14C, pMC, or D14C
        cc.2[,2:3] <- cc.2[,2:3]/1e3
    }
    if(grepl("f", tolower(realm2))) {
        F <- C14toF14C(cc.2[,2], cc.2[,3])
        cc.2[,2:3] <- F
      }
    if(grepl("p", tolower(realm2))) {
      p <- C14topMC(cc.2[,2], cc.2[,3])
      cc.2[,2:3] <- p
    }
    if(grepl("d", tolower(realm2))) {
      F <- C14toF14C(cc.2[,2], cc.2[,3])
      Dmax <- F14CtoD14C(F[,1]+F[,2], t=cc.2[,1])
      D <- F14CtoD14C(F[,1], t=cc.2[,1])
      cc.2[,2:3] <- cbind(D, Dmax-D)
    }

    cc.2 <- cc.2[which(mindat * maxdat == 1),] # limit to the relevant part of the cc only
    if(ka)
      if(grepl("c", tolower(realm2)))
        cc.2 <- cc.2/1e3

    cc2.pol <- cbind(c(cc.2[,cc.cal], rev(cc.2[,cc.cal])), c(cc.2[,2]-cc.2[,3], rev(cc.2[,2]+cc.2[,3])))
  }

  cal.lim <- c(cal1, cal2)
  if(cal.rev)
    cal.lim <- rev(cal.lim)
  if(ka) 
    cal.lim <- cal.lim/1e3
  
  if(!add) { # then prepare plotting parameters
    if(is.na(cal.lab))
      if(ka) {
        if(BCAD) 
          cal.lab <- "ka BC/AD" else
            cal.lab <- "kcal BP"
      } else
        if(BCAD)
          cal.lab <- "BC/AD" else
            cal.lab <- "cal. yr BP"
    if(is.na(c14.lab))
      if(grepl("p", tolower(realm)))
        c14.lab <- "pMC" else
          if(grepl("f", tolower(realm)))
            c14.lab <- expression(F^14*C) else
              if(grepl("d", tolower(realm)))
                c14.lab <- expression(Delta^14*C) else
                  if(ka)
                    c14.lab <- expression(""^14*C~kBP) else
                      c14.lab <- expression(""^14*C~BP)
    if(is.na(c14.lim[1]))
      if(is.na(cc2))
        c14.lim <- range(cc1.pol[,2]) else
          if(add.yaxis)
            c14.lim <- range(cc1.pol[,2]) else
              c14.lim <- range(cc1.pol[,2], cc2.pol[,2])
    if(c14.rev)
      c14.lim <- rev(c14.lim)

    # draw the graph and data
    plot(0, type="n", xlim=cal.lim, xlab=cal.lab, ylim=c14.lim, ylab=c14.lab, bty=bty, ...)
  }

  # add the calibration curve
  polygon(cc1.pol, col=cc1.fill, border=NA) # calibration curve
  lines(cc.1[,cc.cal], cc.1[,2]-cc.1[,3], col=cc1.col)
  lines(cc.1[,cc.cal], cc.1[,2]+cc.1[,3], col=cc1.col)

  # add a second curve?
  if(!is.na(cc2)) {
    if(add.yaxis) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(new=TRUE)
      plot(cc2.pol, type="n", xlim=cal.lim, xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
    }
    polygon(cc2.pol, col=cc2.fill, border=NA) # calibration curve
    lines(cc.2[,cc.cal], cc.2[,2]-cc.2[,3], col=cc2.col)
    lines(cc.2[,cc.cal], cc.2[,2]+cc.2[,3], col=cc2.col)
    if(add.yaxis)
      axis(4, col=cc2.col, col.axis=cc2.col)
    if(length(legend) > 0)
      legend(legend, legend=c(cc1, cc2), text.col=c(cc1.col, cc2.col), bty="n")
  }
  invisible(par("usr")) # for subsequent plot manipulations
}



#' @name calibrate
#' @title Plot individual calibrated dates.
#' @description Calibrate individual 14C dates, plot them and report calibrated ranges.
#' @details
#' Type \code{calibrate()} to see how a date of 2450 +- 50 14C BP gets calibrated (the calibration curve happens to show
#' a plateau around this 14C age). To calibrate a different date, provide its reported mean and error (1
#' standard deviation error as reported by the radiocarbon laboratory) as follows: \code{calibrate(mean, error)},
#' e.g., for a date of 130 +- 10 14C BP, type calibrate\code{(age=130, error=10)} or, shorter, \code{calibrate(130,10)}.
#'
#' In case the date has a reservoir effect or age offset, e.g. of 100 14C years, provide this as follows:
#' \code{calibrate(130, 10, reservoir=100)}. If you want to include an uncertainty for this offset, provide this as follows,
#' e.g., for an uncertainty of 50yr, \code{calibrate(130,10,reservoir=c(100, 50))}.
#' The uncertainty for the age offset will then be added to the error (by taking the square root of the sum
#' of the squared error and the squared offset uncertainty). If the carbon of your sample has mixed marine/terrestrial sources,
#' instead apply the marine offset using mix.curves and calibrate the date using that custom-built curve (cc="mixed").
#'
#' If you prefer to work with, e.g., 68 \% as opposed to the default 95 \% confidence intervals,
#' type: \code{calibrate(130, 10, prob=0.68)} or \code{calibrate(130, 10,, 0.68)} (the commas between the brackets indicate the position of the option;
#' the standard deviation is the fourth option of the \code{calibrate} function). The calibrated distribution can be calculated
#' for every single calendar year (\code{yrsteps=1}) within a wide range of the 14C date. Probabilities below a threshold (default \code{threshold=0.0005}) will be neglected.
#'
#' By default the northern hemisphere terrestrial calibration curve is used (\code{cc=1 or cc1="IntCal20"}).
#' To use alternative curves, use \code{cc=2} (\code{cc2="Marine20"}), \code{cc=3} (\code{cc3="SHCal20C"}),
#' \code{cc=4} (\code{cc4="mixed.14C"}), or specify a postbomb curve (e.g., \code{cc="nh1"}).
#'
#' Calibrate works in cal BP (calendar years before AD 1950) by default, but can work with cal BC/AD through the option \code{BCAD=TRUE}.
#'
#' By default the Gaussian distribution is used to calibrate dates. For use of the t distribution (Christen and Perez 2016) instead,
#' set \code{normal=FALSE} provide values for t.a and t.b (defaults to \code{t.a=3} and \code{t.b=4}).
#'
#' Calibrated distributions are usually reduced to their 68\% or 95\% calibrated ranges, taking into account the asymmetric
#' and multi-peaked shape of these distributions.
#' Calibrated ranges at 68\% will obviously result in narrower confidence intervals, and a perceived higher precision, than 95\% ranges. However, given the often
#' asymmetric and multi-modal nature of calibrated distributions, the probability that the 'true' calendar date
#' lies outside the 1 standard deviation hpd ranges is considerable (c. 32\%). Therefore the use of 95\% calibrated ranges is preferable,
#' and default.
#'
#' Negative radiocarbon ages are calibrated with postbomb curves, but the user needs to tell which curve to use.
#' For example, to use the first of the three northern hemisphere curves, provide the option \code{cc="nh1"}, \code{cc="nh2"}, \code{cc="nh3"},
#' while for southern hemisphere samples, use \code{cc="sh1-2"} or \code{cc="sh3"}.
#'
#' A graph of the calibration is produced, and it can be adapted in several ways.
#' The limits of the horizontal (calendar scale) and vertical (14C scale) axes are calculated automatically
#' but can be changed by providing alternative values for the options \code{cal.lim, C14.lim}.
#' The titles of both axis can be changed by providing alternative titles to \code{cal.lab} and/or \code{C14.lab}. The heights of the distributions of the 14C and calibrated
#' ages can be set to alternative values using \code{dist.height} (default \code{0.3} which plots the distribution up to 30\% of the height of the entire graph).
#' Parameters for white space around the
#' graph can be changed (default \code{mar=c(3.5, 2, 2, 1}) for spacing below, to the left, above and to the right respectively),
#' as can the spacing for the axis labels (\code{mgp=c(2,1,0)}). By default, the axes are connected at the lower left, \code{bty="l"}.
#' Check the R documentation of \code{par()} for more options.
#'
#' The colours of the 14C date, the calibration curve, the distributions, and the highest posterior density (hpd)
#' ranges, can be changed by providing an alternative colour in \code{date.col}, \code{cc.col}, \code{dist.col}, and/or \code{hpd.col}, respectively.
#' The default colours are transparent grey for the dates probability distributions (\code{dist.col=rgb(0,0,0, 0.3)} and \code{sd.col=rgb(0,0,0, 0.5)};
#' change the last value of rgb for different greyscale values), red for the uncalibrated mean and error bars (\code{date.col="red"}),
#' and transparent green for the calibration curve (\code{cc.col=rgb(0, 0.5, 0, 0.7)}). R's rgb() function expects values between \code{0} and \code{1}
#' for red, green and blue, respectively, followed by a value for the semi-transparency (also between 0 and 1). Some graphic devices
#' such as postscript are unable to use transparency; in that case provide different colours or leave the fourth value empty.
#' @param age Mean of the uncalibrated C-14 age.
#' @param error Error of the uncalibrated C-14 age.
#' @param cc Calibration curve for C-14 dates (1, 2, 3, or 4, or, e.g., "IntCal20", "Marine20", "SHCal20", "nh1", "sh3", or "mixed").
#' @param postbomb Whether or not this is a postbomb age. Defaults to FALSE. 
#' @param deltaR Age offset (e.g. for marine samples). Can also be provided as option 'reservoir'.
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation). Can also be provided within option 'reservoir'.
#' @param bombalert Warn if a date is close to the lower limit of the IntCal curve. Defaults to \code{postbomb=TRUE}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). Defaults to c().
#' @param as.F Whether or not to calculate ages in the F14C realm. Defaults to \code{as.F=FALSE}, which uses the C14 realm.
#' @param is.F Use \code{is.F=TRUE} if the date and error are entered as F14C.
#' @param is.pMC Use \code{is.pMC=TRUE} if the date and error are entered as pMC.
#' @param reservoir Reservoir age, or reservoir age and age offset as two values (e.g., \code{reservoir=c(100,10)}).
#' @param prob Probability confidence intervals (between 0 and 1).
#' @param BCAD Use BC/AD or cal BP scale (default cal BP).
#' @param ka Use thousands of years instead of years in the plots and hpd ranges. Defaults to FALSE. 
#' @param draw Whether or not to draw the date. Can be set as FALSE to speed up things
#' @param cal.lab Label of the calendar/horizontal axis. Defaults to the calendar scale, but alternative names can be provided.
#' @param C14.lab Label of the C-14/vertical axis. Defaults to the 14C scale, but alternative names can be provided.
#' @param cal.lim Minimum and maximum of calendar axis (default calculated automatically).
#' @param C14.lim Minimum and maximum of C-14 axis (default calculated automatically).
#' @param cc.col Colour of the lines of the calibration curve. Defaults to semi-transparent dark green; \code{cc.col=rgb(0,.5,0,0.7)}.
#' @param cc.fill Colour of the inner part of the calibration curve. Defaults to semi-transparent dark green; \code{cc.col=rgb(0,.5,0,0.7)}.
#' @param date.col Colour of the "dot-bar" plot of the C14 date. Defaults to \code{date.col="red"}.
#' @param dist.col Colour of the outer lines of the distributions. Defaults to semi-transparent grey, \code{dist.col=rgb(0,0,0,0.2)}.
#' @param dist.fill Colour of the inner part of the distributions. Defaults to semi-transparent grey, \code{dist.col=rgb(0,0,0,0.2)}.
#' @param hpd.fill Colour of the highest posterior density. Defaults to semi-transparent grey, \code{dist.col=rgb(0,0,0,0.3)}.
#' @param dist.height Maximum height of the C14 and calibrated distributions (as proportion of the invisible secondary axes). Defaults to 1.8.
#' @param dist.float The probability distributions float a bit above the axes by default. Can be set to distinct heights of the axes, e.g.: \code{dist.float=c(0.05, 0.1)}, or to \code{dist.float=0}.
#' @param cal.rev Whether or not to reverse the direction of the calendar axis.
#' @param yr.steps Temporal resolution at which C-14 ages are calibrated (in calendar years). By default follows the spacing in the calibration curve.
#' @param cc.resample The IntCal20 curves have different densities (every year between 0 and 5 kcal BP, then every 5 yr up to 15 kcal BP, then every 10 yr up to 25 kcal BP, and then every 20 yr up to 55 kcal BP). If calibrated ages span these density ranges, their drawn heights can differ, as can their total areas (which should ideally all sum to the same size). To account for this, resample to a constant time-span, using, e.g., cc.resample=5 for 5-yr timespans.
#' @param threshold Below which value should probabilities be excluded from calculations.
#' @param edge How to treat dates are at or beyond the edge of the calibration curve. If dates are truncated, a warning is given. If they lie beyond the calibration curve, an error is given.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2016).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param rounded Rounding of the percentages of the reported hpd ranges. Defaults to 1 decimal.
#' @param every Yearly precision of the hpd ranges (defaults to \code{every=1}).
#' @param extend.range Range by which the axes are extended beyond the data limits. Defaults to 5\%.
#' @param legend.cex Size of the font of the legends. Defaults to 0.8.
#' @param legend1.loc Where the first legend (with the calibration curve name and the uncalibrated date) is plotted. Defaults to topleft.
#' @param legend2.loc Where the second legend (with the hpd ranges) is plotted. Defaults to topright.
#' @param print.truncate.warning Whether or not a truncation warning is printed on the plot. Defaults to \code{print.truncate.warning=TRUE}.
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted).
#' @param mar Plot margins (amount of white space along edges of axes 1-4).
#' @param xaxs Whether or not to extend the limits of the horizontal axis. Defaults to \code{xaxs="i"} which does not extend the limits.
#' @param yaxs Whether or not to extend the limits of the vertical axis. Defaults to \code{yaxs="i"} which does not extend the limits.
#' @param bty Draw a box around the graph ("n" for none, and "l", "7", "c", "u", "]" or "o" for correspondingly shaped boxes).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param cc.er The error of the calibration curve. Only used for plotting the uncalibrated C14 distribution, which by default only shows the date's uncertainty (the calibration curve uncertainty is indeed taken into account during calibration). If known, the calibration curve's error can be added.
#' @param every Yearly precision (defaults to \code{every=1}).
#' @param ... Other plotting parameters.
#' @return A graph of the raw and calibrated C-14 date, the calibrated ranges and, invisibly, the calibrated distribution and hpd ranges.
#' @examples
#' calibrate()
#' calibrate(130, 10)
#' cal <- calibrate(2550, 20, reservoir=100)
#' cal; plot(cal[[1]])
#' calibrate(130, 10, prob=0.68)
#' calibrate(age=130, error=10, BCAD=TRUE)
#' calibrate(4450, 40, reservoir=c(100, 50))
#' @export
calibrate <- function(age=2450, error=50, cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, bombalert=TRUE, thiscurve=c(), as.F=FALSE, is.F=FALSE, is.pMC=FALSE, reservoir=0, prob=0.95, BCAD=FALSE, ka=FALSE, draw=TRUE, cal.lab=c(), C14.lab=c(), cal.lim=c(), C14.lim=c(), cc.col=rgb(0,.5,0,0.7), cc.fill=rgb(0,.5,0,0.7), date.col="red", dist.col=rgb(0,0,0,0.3), dist.fill=rgb(0,0,0,0.3), hpd.fill=rgb(0,0,0,0.3), dist.height=0.3, dist.float=c(.01, .01), cal.rev=FALSE, yr.steps=FALSE, cc.resample=5, threshold=0.0005, edge=TRUE, normal=TRUE, t.a=3, t.b=4, rounded=1, every=1, extend.range=.05, legend.cex=0.8, legend1.loc="topleft", legend2.loc="topright", print.truncate.warning=TRUE, mgp=c(2,1,0), mar=c(3,3,1,1), xaxs="i", yaxs="i", bty="l", cc.dir=NULL, cc.er=0, ...) {
  
  if(is.F && is.pMC)
    stop("Cannot have both is.F=TRUE and is.PMC=TRUE.")
  
  age <- age - deltaR
  error <- sqrt(error^2 + deltaSTD^2)
  
  # read the data
  age <- age-reservoir[1]
  if(length(reservoir) > 1)
    error <- sqrt(error^2 + reservoir[2]^2)
   
  # check if the date is covered by the curve
  beyond <- FALSE
  youngest.cc <- c(95,603,118,0,0) # youngest C14 ages of, respectively, IntCal20, Marine20, SHCal20, and extra entries
  if(any(cc %in% 1:4)) {
    if(is.F) {
      beyond <- (age + (3*error)) > C14toF14C(youngest.cc)[cc]
    } else 
        if(is.pMC) {
          beyond <- (age + (3*error)) > C14topMC(youngest.cc)[cc]
        } else
      beyond <- (age - (3*error)) < youngest.cc[cc]
  }

  if(bombalert) {
    if(beyond) { # at or beyond younger IntCal limit
      if(!postbomb) # note that there are no postbomb curves for Marine20
        if(!(cc %in% c("nh1", "nh2", "nh3", "sh1-2", "sh3", "nh1_monthly", "nh2_monthly", "nh3_monthly", "sh1-2_monthly", "sh3_monthly")))
          stop("This appears to be a postbomb age (or is close to being one). Please provide a postbomb curve")
    Cc <- rintcal::glue.ccurves(cc, postbomb, cc.dir) # doesn't do resample
  } else {
      if(postbomb > 0) # postbomb has been defined
        Cc <- rintcal::glue.ccurves(cc, postbomb, cc.dir) else
          Cc <- rintcal::ccurve(cc, postbomb=FALSE, cc.dir)
    }
  } else
    Cc <- rintcal::ccurve(cc, postbomb=postbomb, cc.dir)
  cc.cal <- 1
  if(BCAD) {
    Cc[,4] <- calBPtoBCAD(Cc[,1])
    cc.cal <- 4
  }
  if(length(thiscurve) > 0) # then forget the above
    Cc <- thiscurve

  if(is.F) # then also put Cc on the F scale
    Cc <- cbind(Cc[,1], C14toF14C(Cc[,2], Cc[,3]), Cc[,cc.cal])
  if(is.pMC) # then also put Cc on the F scale
    Cc <- cbind(Cc[,1], C14topMC(Cc[,2], Cc[,3]), Cc[,cc.cal])
  
  # warn/stop if the date lies (partly) beyond the calibration curve
  if(is.pMC)
    asF <- c(age, error)/100 else
      if(is.F)
        asF <- c(age, error) else
          asF <- C14toF14C(age, error)
    # maxcc <- rintcal::ccurve(cc, postbomb=postbomb, cc.dir, as.F=TRUE)
    # maxcc <- maxcc[nrow(maxcc),2:3]
    if(asF[1] < (2*asF[2])) { # according to background conventions
    msg <- paste("age < (2*er) in F14C space; this date is at background and should NOT be calibrated!")
    if(edge)
      stop(msg) else
        message(msg)
  }
  truncate.warning <- FALSE

  border <- 0
  Age <- age; Error <- error
  if(Age-2*Error < min(Cc[,2]-2*Cc[,3]))
    if(Age+2*Error > min(Cc[,2]+2*Cc[,3]))
      border <- 1 else border <- 2
  if(Age+2*Error > max(Cc[,2]-2*Cc[,3]))
    if(Age-2*Error < min(Cc[,2]+2*Cc[,3]))
      border <- 1 else border <- 2

  if(border == 1)
    if(edge)
      message("Date falls partly beyond calibration curve and will be truncated!")
  if(border == 2)
    if(edge)
      stop("Cannot calibrate dates beyond the calibration curve!")
  if(border > 0)
    truncate.warning <- TRUE

  C14.dist <- caldist(age, sqrt(error^2 + cc.er^2), cc=0, BCAD=FALSE, postbomb=FALSE) # just to draw a normal dist
  cal.dist <- caldist(age, error, cc=cc, yrsteps=yr.steps, threshold=threshold,
  normal=normal, is.F=is.F, is.pMC=is.pMC, as.F=as.F, t.a=t.a, t.b=t.b, postbomb=postbomb, cc.dir=cc.dir, BCAD=BCAD)

  # copy entries at edges of calibrated hpds, to ensure rectangular polygons
  if(draw) {

    # calculate limits
    if(length(cal.lim) == 0) {
      cal.lim <- range(cal.dist[,1])
    lims <- cal.lim
    cal.lim <- rev(extendrange(cal.lim, f=extend.range))
    if(BCAD)
      cal.lim <- rev(cal.lim)
    if(cal.rev)
      cal.lim <- rev(cal.lim)
    }  
    if(length(C14.lim) == 0) {
      if(BCAD) {
        cc.min <- max(1, min(which(Cc[,cc.cal] <= max(cal.lim))))
        cc.max <- min(nrow(Cc), max(which(Cc[,cc.cal] >= min(cal.lim))))
      } else {  
          cc.min <- max(1, min(which(Cc[,1] >= min(cal.lim))))
          cc.max <- min(nrow(Cc), max(which(Cc[,1] <= max(cal.lim))))
        }
      # we don't need the entire calibration curve
      Cc <- Cc[cc.min:cc.max,]
      if(BCAD)
        cc.lim <- extendrange(c(Cc[,2]-Cc[,3], Cc[,2]+Cc[,3]), f=extend.range) else
          cc.lim <- extendrange(c(Cc[,2]-Cc[,3], Cc[,2]+Cc[,3], C14.dist[,1]), f=extend.range)
    } else {
       cc.min <- max(1, which(Cc[,2] >= min(C14.lim)))
       cc.max <- min(nrow(Cc), which(Cc[,2] <= max(C14.lim)))
       Cc <- Cc[cc.min:cc.max,]
       cc.lim <- range(C14.lim)
    }
    ccpol <- cbind(c(Cc[,cc.cal], rev(Cc[,cc.cal])), c(Cc[,2]-Cc[,3], rev(Cc[,2]+Cc[,3])))

    if(length(dist.float) == 1)
      dist.float[2] <- dist.float[1] 
    callim <- cal.lim[1]+dist.float[2]*(cal.lim[2]-cal.lim[1])
    cclim <- cc.lim[1]+dist.float[1]*(cc.lim[2]-cc.lim[1])
  
    # adapt axis titles, labels and hpds if BCAD and/or ka
    if(length(cal.lab) == 0)
      if(ka) 
        cal.lab <- ifelse(BCAD, "k BC/AD", "kcal BP") else 
          cal.lab <- ifelse(BCAD, "BC/AD", "cal BP")
    if(length(C14.lab) == 0)
      if(is.F)
        C14.lab <- expression(F^14*C) else
          if(is.pMC)
            C14.lab <- "pMC" else
              if(ka)
                C14.lab <- expression(""^14*C~kBP) else
                  C14.lab <- expression(""^14*C~BP)
    xaxt <- ifelse(BCAD || ka, "n", "s")
    yaxt <- ifelse(ka, "n", "s")

    plot(0, type="n", xlim=cal.lim, ylim=cc.lim, xlab=cal.lab, ylab=C14.lab, xaxt="n", yaxt="n", xaxs=xaxs, yaxs=yaxs, bty=bty, mgp=mgp, mar=mar)
    if(ka) {
      axis(1, pretty(cal.lim), labels=pretty(cal.lim/1e3))
      axis(2, pretty(cc.lim), labels=pretty(cc.lim/1e3))
    } else {
       axis(1)
       axis(2)
    }

    # draw the data
    coors <- par('usr')
    polygon(ccpol, border=cc.col, col=cc.fill)
    draw.dist(C14.dist, on.y=TRUE, x.pos=callim, as.unit=FALSE, fraction=dist.height, mirror=FALSE, up=TRUE, peak=TRUE, prob=prob, BCAD=BCAD, hpd.border=NA, hpd.col=hpd.fill, dist.col=dist.fill, dist.border=dist.col)
    dot <- ifelse(cal.rev, min(lims), callim)
    points(dot, age, col=date.col, pch=20)
    segments(dot, age-error, dot, age+error, col=date.col)
    hpds <- draw.dist(cal.dist, on.y=FALSE, y.pos=cclim, as.unit=FALSE, fraction=dist.height, mirror=FALSE, up=TRUE, peak=TRUE, prob=prob, BCAD=BCAD, hpd.border=NA, hpd.col=hpd.fill, dist.col=dist.fill, dist.border=dist.col)

    # legends
    if(cc == 1)
      cc <- "IntCal20 " else
        if(cc == 2)
          cc <- "Marine20 " else
            if(cc == 3)
              cc <- "SHCal20 "
    legend(legend1.loc, legend=c(cc, paste(age, "\u00B1", error)), text.col=c(cc.col, 1),  ncol=1, bty="n", cex=legend.cex)
    legend(legend2.loc, legend=rbind(c("from", "to", "%"), cbind(round(hpds, rounded))), ncol=3, bty="n", cex=legend.cex)
  }
  
  if(truncate.warning)
    if(print.truncate.warning)
      legend("right", "Warning! Date truncated ", text.col="red", bty="n", cex=legend.cex, xjust=1)
  
  invisible(list(calib=cal.dist, hpd=hpds))
}



# internal function to draw distributions
draw.dist <- function(dist, on.y=FALSE, rotate.axes=FALSE, mirror=FALSE, up=TRUE, hpd=TRUE, ka=FALSE, prob=0.95, prob.round=1, BCAD=FALSE, x.pos=c(), y.pos=c(), ex=1, peak=TRUE, as.unit=FALSE, fraction=0.1, dist.col=rgb(0,0,1,0.3), dist.border=rgb(0,0,1,0.3), hpd.col=rgb(0,0,1,0.3), hpd.border=rgb(0,0,1,0.3)) {

  ages0 <- c(dist[1,1], dist[,1], dist[nrow(dist),1])
  agesmirror <- c(dist[,1], rev(dist[,1]))
  max.peak <- max(dist[,2])
  if(peak) # then set peak height to 1
    dist[,2] <- dist[,2] / max.peak
  dist0 <- ex*c(0, dist[,2], 0)
  distmirror <- ex/2*c(dist[,2], -1*rev(dist[,2]))

  # which direction to draw the distribution(s)
  xy <- par('usr')
  if(on.y) {
	leftright <- fraction*(xy[2] - xy[1]) # from left to right corner
    if(as.unit)
      leftright <- ifelse(leftright >= 0, 1, -1)
    add <- ifelse(up, 1, -1)*ex*leftright
  } else {
    downup <- fraction*(xy[4] - xy[3]) # from bottom to top corner
    if(as.unit)
      downup <- ifelse(downup >= 0, 1, -1)
    add <- ifelse(up, 1, -1)*ex*downup 
  }

  if(on.y) { # draw on the left axis
    if(length(x.pos) == 0)
      x.pos <- xy[1] # left extreme of x axis
    if(mirror)
      pol <- cbind(x.pos+distmirror, agesmirror) else # x.pos was y.pos
        pol <- cbind(x.pos+add*dist0, ages0)
  } else {
      if(length(y.pos) == 0)
        y.pos <- xy[3] # bottom extreme of y axis
      if(mirror)
        pol <- cbind(agesmirror, y.pos+distmirror) else
          pol <- cbind(ages0, y.pos+add*dist0)
    }
    polygon(pol, col=dist.col, border=dist.border)

  if(hpd) {
    hpds <- hpd(dist, return.raw=TRUE, BCAD=BCAD, prob=prob, prob.round=prob.round, ka=ka)
    hpds.raw <- hpds$calib
    hpds <- hpds$hpds

    for(i in 1:nrow(hpds)) {
      thisdist <- (dist[,1] >= min(hpds[i,1:2])) * (dist[,1] <= max(hpds[i,1:2]))
      thisdist <- dist[which(thisdist == 1),]

      if(length(thisdist) == 2) # then just 1 row
        thisdist <- rbind(c(thisdist[1],0), thisdist, c(thisdist[1], 0))
      ages0 <- c(thisdist[1,1], thisdist[,1], thisdist[nrow(thisdist),1])
      agesmirror <- c(thisdist[,1], rev(thisdist[,1]))
      dist0 <- ex*c(0, thisdist[,2], 0)
      distmirror <- ex/2*c(thisdist[,2], -1*rev(thisdist[,2]))

      if(on.y) {
        if(length(x.pos) == 0)
          x.pos <- par('usr')[1] # left extreme of x axis
      if(mirror)
        pol <- cbind(y.pos+distmirror, agesmirror) else
          pol <- cbind(x.pos+add*dist0, ages0) 
      } else {
        if(length(y.pos) == 0)
          y.pos <- par('usr')[3] # bottom extreme of y axis
        if(mirror)
          pol <- cbind(agesmirror, y.pos+distmirror) else
            pol <- cbind(ages0, y.pos+add*dist0)
      }
    polygon(pol, col=dist.col, border=dist.border)
    }
  invisible(hpds)
  }
}



#' @name draw.dates
#' @title add calibrated distributions to a plot.
#' @description Add individual or multiple calibrated dates to a plot.
#' @return A plot of the (calibrated) dates
#' @param age Mean of the uncalibrated C-14 age (or multiple ages).
#' @param error Error of the uncalibrated C-14 age (or ages).
#' @param depth Depth(s) of the date(s). Defaults to their relative positions if no depths are provided.
#' @param cc Calibration curve for C-14 dates (1, 2, 3, or 4, or, e.g., "IntCal20", "Marine20", "SHCal20", "nh1", "sh3", or "mixed"). If there are multiple dates but all use the same calibration curve, one value can be provided. 
#' @param postbomb Whether or not this is a postbomb age. Defaults to FALSE. 
#' @param deltaR Age offset (e.g. for marine samples). Can also be provided as option 'reservoir'.
#' @param deltaSTD Uncertainty of the age offset (1 standard deviation). Can also be provided within option 'reservoir'.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error). Defaults to c().
#' @param oncurve Whether or not to plot the calibration curve and then plot the dates onto this curve. Defaults to FALSE.
#' @param realm If oncurve is used, by default the calibration curve is plotted in the C14 age realm. Alternatively, it can be provided as \code{realm="F14C"} or \code{realm="pMC"}  
#' @param reservoir Reservoir age, or reservoir age and age offset.
#' @param normal Use the normal distribution to calibrate dates (default TRUE). The alternative is to use the t model (Christen and Perez 2009).
#' @param t.a Value a of the t distribution (defaults to 3).
#' @param t.b Value b of the t distribution (defaults to 4).
#' @param prob Probability confidence intervals (between 0 and 1).
#' @param threshold Report only values above a threshold. Defaults to \code{threshold=0.001}.
#' @param BCAD Use BC/AD or cal BP scale (default cal BP).
#' @param draw.hpd Whether or not to draw the hpd ranges as a line
#' @param hpd.col Colour of the hpd rectangle for all dates or radiocarbon dates
#' @param hpd.border Colour of the border of the hpd intervals. Not drawn by default.
#' @param cal.hpd.col Colour of the hpd rectangle for cal BP dates
#' @param rounded Rounding for probabilities of reported hpd ranges. Defaults to 1 decimal.
#' @param every Yearly precision of hpds (defaults to \code{every=1}).
#' @param mirror Plot distributions mirrored, a bit like a swan. Confuses some people but looks nice to the author so is the default.
#' @param up If mirror is set to FALSE, the distribution can be plotted facing upwards or downwards.
#' @param draw.base By default, the base of the calibrated distributions is plotted. This can be avoided by supplying \code{draw.base=FALSE} as an option.
#' @param col Colour of the inside of the distribution
#' @param border Colour of the border of the distribution
#' @param cal.col Colour of the inside of distribution of non-radiocarbon dates that didn't need calibration
#' @param cal.border Colour of the border of the distribution of non-radiocarbon dates that didn't need calibration
#' @param add Whether or not to add the dates to an existing plot. If set to FALSE (default), a plot will be set up.
#' @param ka Whether or not to plot ages as thousands of years. Defaults to \code{ka=FALSE}.
#' @param rotate.axes By default, the calendar age axis is plotted on the horizontal axis, and depth/position on the vertical one. Use \code{rotate.axes=TRUE} to rotate the axes.
#' @param ex Exaggeration of the height of the distribution, defaults to \code{ex=1}.
#' @param normalise If TRUE, the age distributions are normalised by plotting each distribution with the same total area. Precise dates will therefore peak higher than less precise dates (default). If \code{normalise=FALSE}, the peak of each date will be drawn at the same height.
#' @param cc.col Colour of the calibration curve. Default semi-transparent darkgreen.
#' @param cc.border Colour of the edges of the calibration curve. Default semi-transparent darkgreen.
#' @param cc.resample The IntCal20 curves have different densities (every year between 0 and 5 kcal BP, then every 5 yr up to 15 kcal BP, then every 10 yr up to 25 kcal BP, and then every 20 yr up to 55 kcal BP). If calibrated ages span these density ranges, their drawn heights can differ, as can their total areas (which should ideally all sum to the same size). To account for this, resample to a constant time-span, using, e.g., \code{cc.resample=5} for 5-yr timespans.
#' @param age.lab Title of the calendar axis (if present)
#' @param age.lim Limits of the calendar axis (if present)
#' @param age.rev Reverse the age axis. Defaults to TRUE
#' @param d.lab Title of the vertical axis (if present)
#' @param d.lim Limits of the vertical axis (if present)
#' @param d.rev Reverse the y-axis. Defaults to TRUE
#' @param labels Add labels to the dates. Empty by default.
#' @param label.x Horizontal position of the date labels. By default draws them before the youngest age (1), but can also draw them after the oldest age (2), or above its mean (3). 
#' @param label.y Vertical positions of the depths/labels. Defaults to 0 (or 1 if label.x is 3 or 4).
#' @param label.offset Offsets of the positions of the depths/labels, giving the x and y offsets. Defaults to c(0,0).
#' @param label.cex Size of labels. 
#' @param label.col Colour of the labels. Defaults to the colour given to the borders of the dates.
#' @param label.adj  Justification of the labels. Follows R's adj option: A value of "0" produces left-justified text, "0.5" (the default) centered text and "1" right-justified text.
#' @param label.rot Rotation of the label. 0 by default (horizontal).
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param dist.res Resolution of the distribution polygons. Defaults to \code{dist.res=100}.
#' @param ... Additional plotting options
#' @examples
#'   plot(0, xlim=c(500,0), ylim=c(0, 2))
#'   draw.dates(130, 20, depth=1) 
#'   x <- sort(runif(10, 1000, 10000)) # draw 10 random calendar ages
#'   cc <- rintcal::ccurve() # get the calibration curve
#'   y <- approx(cc[,1], cc[,2], x)$y # find the IntCal 14C ages
#'   er <- .01 * y
#'   draw.dates(y, er, 1:length(x))
#'   # or draw on the calibration curve
#'   draw.dates(y, er, y, d.lab="Radiocarbon age (BP)")
#'   draw.ccurve(add=TRUE, cc1.col=rgb(0,.5,0,.5))
#' @export
draw.dates <- function(age, error, depth=c(), cc=1, postbomb=FALSE, deltaR=0, deltaSTD=0, thiscurve=c(), oncurve=FALSE, realm="C", reservoir=c(), normal=TRUE, t.a=3, t.b=4, prob=0.95, threshold=.001, BCAD=FALSE, draw.hpd=TRUE, hpd.border=NA, hpd.col=rgb(0,0,1,.7), cal.hpd.col=rgb(0, 0.5, 0.5, 0.35), rounded=0.1, every=1, mirror=TRUE, up=TRUE, draw.base=TRUE, col=rgb(0,0,1,.3), border=rgb(0,0,1,.5), cal.col=rgb(0, 0.5, 0.5, 0.35), cal.border=rgb(0, 0.5, 0.5, 0.35), add=FALSE, ka=FALSE, rotate.axes=FALSE, ex=.8, normalise=TRUE, cc.col=rgb(0,.5,0,.5), cc.border=rgb(0,.5,0,.5), cc.resample=5, age.lab=c(), age.lim=c(), age.rev=FALSE, d.lab=c(), d.lim=c(), d.rev=TRUE, labels=c(), label.x=1, label.y=c(), label.cex=0.8, label.col=border, label.offset=c(0,0), label.adj=c(1,0), label.rot=0, cc.dir=NULL, dist.res=100, ...) {

  age <- age - deltaR
  error <- sqrt(error^2 + deltaSTD^2)

  if(length(reservoir) > 0) {
    age <- age - reservoir[1]
    if(length(reservoir) > 1)
      error <- sqrt(error^2 + reservoir[2]^2)
  }
  if(length(depth) == 0)
    depth <- 1:length(age)    

  # deal with multiple dates
  if(length(age) > 1) {
    if(length(cc) == 1)
      cc <- rep(cc, length(age)) 
    if(length(postbomb) == 1)
      postbomb <- rep(postbomb, length(age))
    if(length(hpd.col) == 1)
      hpd.col <- rep(hpd.col, length(age))
    if(length(hpd.border) == 1)
      hpd.border <- rep(hpd.border, length(age))
    if(length(col) == 1)
      col <- rep(col, length(age))
    if(length(border) == 1)
      border <- rep(border, length(age))
    if(0 %in% cc) { # then we've got cal BP dates
      these <- which(cc==0)
      hpd.col[these] <- cal.hpd.col
      hpd.border[these] <- cal.hpd.col
      col[these] <- cal.col
      border[these] <- cal.border
    }
  }

  ages <- array(NA, dim=c(dist.res, length(age))) # later fill with years
  probs <- ages # later fill with probs
  mx <- rep(0, length(age))
  hpds <- list()
  for(i in 1:length(age)) {
    tmp <- caldist(age[i], error[i], cc=cc[i], postbomb=postbomb[i], normal=normal, t.a=t.a, t.b=t.b, normalise=normalise, thiscurve=thiscurve, cc.resample=cc.resample, threshold=threshold, BCAD=BCAD, cc.dir=cc.dir)

    tmp <- approx(tmp[,1], tmp[,2], seq(min(tmp[,1]), max(tmp[,1]), length=dist.res))
    if(normalise)
      tmp <- cbind(tmp$x, tmp$y/sum(tmp$y)) else
        tmp <- cbind(tmp$x, tmp$y/max(tmp$y))

    mx[i] <- max(tmp[,2])
    ages[,i] <- tmp[,1]
    probs[,i] <- tmp[,2]
  }

  if(normalise)
    probs <- probs/max(probs) # otherwise they are very low
  if(ka)
    ages <- ages/1e3

  ages <- cbind(ages)
  probs <- cbind(probs)
  if(oncurve) {
    if(ka)
      depth <- age/1e3 else
        depth <- age  
    if(grepl("p", tolower(realm)))
      depth <- C14topMC(depth)
    if(grepl("f", tolower(realm)))
      depth <- C14toF14C(depth)
  }

  if(!add) {
    if(length(age.lab) == 0)
      if(ka) {
        age.lab <- ifelse(BCAD, "k BC/AD", "kcal BP")
      } else
        age.lab <- ifelse(BCAD, "BC/AD", "cal BP")

    if(length(d.lab) == 0)
      if(oncurve) {
        if(grepl("p", tolower(realm)))
          d.lab <- "pMC" else 
         if(grepl("f", tolower(realm)))
           d.lab <- "F14C" else
           if(ka)
             d.lab <- "C14 kBP" else
               d.lab <- "C14 BP" 
      } else
          d.lab <- "depth"

    if(length(age.lim) == 0)
      age.lim <- extendrange(ages)
    if(!BCAD)
      age.lim <- age.lim[2:1]
    if(length(d.lim) == 0)
      d.lim <- extendrange(c(depth-ex, depth+ex))
    if(d.rev)
      d.lim <- rev(d.lim)
    if(age.rev)
      age.lim <- rev(age.lim)
    if(oncurve)
      age.lim <- rev(age.lim)
    if(rotate.axes)
      plot(0, type="n", ylim=age.lim, ylab=age.lab, xlim=d.lim, xlab=d.lab, ...) else
        plot(0, type="n", xlim=age.lim, xlab=age.lab, ylim=d.lim, ylab=d.lab, ...)
    if(oncurve) {
      if(length(thiscurve) > 0)
        cc <- thiscurve else
          cc <- rintcal::ccurve(cc=cc[1], postbomb=postbomb[1])
        if(BCAD)
          cc[,1] <- calBPtoBCAD(cc[,1])
        if(grepl("p", tolower(realm)))
          cc[,2:3] <- C14topMC(cc[,2], cc[,3])
        if(grepl("f", tolower(realm)))
          cc[,2:3] <- C14toF14C(cc[,2], cc[,3])
        if(ka)
          cc <- cc/1e3
      ccpol <- cbind(c(cc[,1], rev(cc[,1])), c(cc[,2]-cc[,3], rev(cc[,2]+cc[,3])))
      polygon(ccpol, col=cc.col, border=cc.border)
    }
  } else {
      ax.lm <- par("usr")
      if(rotate.axes)
        d.lim <- ax.lm[1:2] else
          d.lim <- ax.lm[3:4]
  }

#   scaling <- c() # to plot multiple dates without overlap
#   if(length(age) > 1) {
#     alldepths <- sort(c(d.lim, depth))
#     if(mirror)
#       for(i in 2:(length(alldepths)-1)) {
#         available <- min(abs(diff(alldepths[(i-1):(i+1)])))
#         scaling <- c(scaling, ex * available / mx[i])
#       } else
#           for(i in 2:length(alldepths)) {
#             available <- min(abs(diff(alldepths[(i-1):(i)])))
#             scaling <- c(scaling, ex * available / mx[i])
#           }
#   } else
#     scaling <- 1

  # now draw the dates
  for(i in 1:length(age)) {
    if(draw.base) {
      if(rotate.axes)
        draw.dist(cbind(ages[,i], probs[,i]), up=up, mirror=mirror, ka=ka, on.y=TRUE, rotate.axes=TRUE, ex=ex, x.pos=depth[i], hpd=draw.hpd, prob=prob, hpd.col=hpd.col[i], hpd.border=hpd.border[i], dist.col=col, dist.border=border) else
          draw.dist(cbind(ages[,i], probs[,i]), up=up, mirror=mirror, ka=ka, on.y=FALSE, ex=ex, y.pos=depth[i], hpd=draw.hpd, prob=prob, hpd.col=hpd.col[i], hpd.border=hpd.border[i], dist.col=col, dist.border=border)
    } else
        if(rotate.axes)
          lines(depth[i]+ex*probs[,i], ages[,i], col=col[i]) else
            lines(ages[,i], depth[i]+ex*probs[,i], col=col[i])

    if(length(labels) > 0) {
      xx <- c(min(ages[,i]), max(ages[,i]), mean(ages[,i]))
      if(!BCAD) xx <- xx[c(2,1,3)]
      x <- xx[label.x]
      if(length(label.y) == 0) {
        y <- depth[i]
        if(label.x > 2)
          ifelse(up, y <- y+1, y <- y-1)
      }

    if(rotate.axes)
       text(y+label.offset[1], x+label.offset[2], labels[i], cex=label.cex, col=label.col, adj=label.adj, srt=label.rot) else
         text(x+label.offset[2], y+label.offset[1], labels[i], cex=label.cex, col=label.col, adj=label.adj, srt=label.rot)
    }
  }
  invisible(list(ages=ages, probs=probs))
}



#' @name draw.D14C
#' @title Draw d14C and the calibration curve.
#' @description Draw a proxy of the atmospheric 14C concentration (d14C) as well as the calibration curve.
#' @return A plot of d14C and the calibration curve
#' @param cal1 First calendar year for the plot. Defaults to youngest calendar age of the calibration curve
#' @param cal2 Last calendar year for the plot. Defaults to oldest calendar age of the calibration curve
#' @param cc The calibration curve to use. Defaults to IntCal20
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param mar Plot margins (amount of white space along edges of axes 1-4).
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted).
#' @param xaxs Whether or not to extend the limits of the horizontal axis. Defaults to \code{xaxs="r"} which extends it by R's default.
#' @param yaxs Whether or not to extend the limits of the vertical axis. Defaults to \code{yaxs="r"} which extends it by R's default.
#' @param bty Draw a box around the graph ("n" for none, and "l", "7", "c", "u", "]" or "o" for correspondingly shaped boxes).
#' @param ka Use kcal BP (and C14 kBP). Defaults to FALSE.
#' @param cal.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}), or to \code{age.lab="kcal BP"} etc. if ka=TRUE.
#' @param cal.rev Reverse the calendar axis (defaults to FALSE).
#' @param C14.lab Label for the C-14 axis. Defaults to 14C BP (or 14C kBP if ka=TRUE).
#' @param C14.lim Limits for the C-14 axis. Calculated automatically by default.
#' @param cc.col Colour of the calibration curve (fill).
#' @param cc.border Colour of the calibration curve (border).
#' @param D14C.lab Label for the D14C axis.
#' @param D14C.lim Axis limits for the D14C axis. Calculated automatically by default.
#' @param D14C.col Colour of the D14C curve (fill).
#' @param D14C.border Colour of the D14C curve (border).
#' @examples
#'   draw.D14C()
#'   draw.D14C(30e3, 55e3, ka=TRUE)
#'   draw.D14C(cc=rintcal::ccurve("NH1_monthly"), BCAD=TRUE)
#' @export
draw.D14C <- function(cal1=c(), cal2=c(), cc=rintcal::ccurve(), BCAD=FALSE, mar=c(4,4,1,4), mgp=c(2.5,1,0), xaxs="r", yaxs="r", bty="u", ka=FALSE, cal.lab=c(), cal.rev=FALSE, C14.lab=c(), C14.lim=c(), cc.col=rgb(0,.5,0,.5), cc.border=rgb(0,.5,0,.5), D14C.lab=c(), D14C.lim=c(), D14C.col=rgb(0,0,1,.5), D14C.border=rgb(0,0,1,.5)) {
  cc.cal <- 1
  if(BCAD) {
    cc[,4] <- 1950 - cc[,1] # add a column
    cc.cal <- 4
  if(length(cal1) == 0)
    cal1 <- min(cc[,cc.cal]) else
      cal1 <- cc[min(1, which(cc[,cc.cal] >= cal1)),1] # find the cal BP values
  if(length(cal2) == 0)
    cal2 <- min(cc[,cc.cal]) else
      cal2 <- cc[min(nrow(cc), which(cc[,cc.cal] <= cal1)),1]
  } else { 
      if(length(cal1) == 0)
        cal1 <- min(cc[,cc.cal])
      if(length(cal2) == 0)
        cal2 <- max(cc[,cc.cal])
    }

  yrmin <- max(1, which(cc[,1] <= min(cal1, cal2)))
  yrmax <- min(nrow(cc), which(cc[,1] >= max(cal1, cal2)))
  cc <- cc[yrmin:yrmax,]
  cc.Fmin <- C14toF14C(cc[,2]+cc[,3])
  cc.Fmax <- C14toF14C(cc[,2]-cc[,3])
  cc.D14Cmin <- F14CtoD14C(cc.Fmin, t=cc[,1])
  cc.D14Cmax <- F14CtoD14C(cc.Fmax, t=cc[,1])
  op <- par(mar=mar, bty=bty, mgp=mgp, xaxs=xaxs, yaxs=yaxs)
  on.exit(par(op))
  kyr <- ifelse(ka, 1e3, 1)
  cal.lim <- range(cal1, cal2)/kyr
  if(BCAD)
    cal.lim <- 1950/kyr - cal.lim 
  if(cal.rev)
    cal.lim <- rev(cal.lim)
  if(length(cal.lab) == 0)
    if(BCAD)
      cal.lab <- ifelse(kyr > 1, "kcal BC/AD", "cal BC/AD") else
        cal.lab <- ifelse(kyr > 1, "kcal BP", "cal BP")
  if(length(D14C.lab) == 0)
    D14C.lab <- expression(Delta^14*C*" (\u{2030})") # assuming UTF8 capabilities
  if(length(C14.lab) == 0)
    C14.lab <- ifelse(kyr > 1, expression(""^14*C~kBP), expression(""^14*C~BP))
  if(length(D14C.lim) == 0)
    D14C.lim <- range(cc.D14Cmin, cc.D14Cmax)
  if(length(C14.lim) == 0)
    C14.lim <- range((cc[,2]-cc[,3]), (cc[,2]+cc[,3]))
  if(ka)
    C14.lim <- C14.lim/1e3
  plot(cc[,cc.cal]/kyr, cc.D14Cmax, type="n", xlab=cal.lab, yaxt="n", ylab="", xlim=cal.lim, ylim=D14C.lim)
  pol.D14C <- cbind(c(cc[,cc.cal]/kyr, rev(cc[,cc.cal]/kyr)), c(cc.D14Cmin, rev(cc.D14Cmax)))
  polygon(pol.D14C, col=D14C.col, border=D14C.border)
  axis(2, col=D14C.border, col.axis=D14C.border)
  mtext(D14C.lab, 2, mgp[1], col=D14C.border)

  op <- par(new=TRUE)
  on.exit(par(op))
  plot(cc[,cc.cal]/kyr, (cc[,2]+cc[,3])/kyr, type="n", xaxt="n", yaxt="n", col=4, xlim=cal.lim, xlab="", ylab="", ylim=C14.lim)
  pol.cc <- cbind(c(cc[,cc.cal]/kyr, rev(cc[,cc.cal]/kyr)), c((cc[,2]+cc[,3])/kyr, rev((cc[,2]-cc[,3])/kyr)))
  polygon(pol.cc, col=cc.col, border=cc.border)

  axis(4, col=cc.border, col.axis=cc.border)
  mtext(C14.lab, 4, mgp[1], col=cc.border)
}



#' @name draw.contamination
#' @title Draw contamination impacts
#' @description Show how contamination with different fractions of modern carbon affect observed C-14 ages.
#' @return A plot of real and observed (contamination-impacted) C14 ages.
#' @param from Minimum 14C age for the plot. Defaults to 0
#' @param to Maximum 14C age for the plot. Defaults to 50e3.
#' @param ka Use C14 kBP. Defaults to TRUE.
#' @param age.res Resolution of age scale. Defaults to 500, which results in smooth curves. Higher numbers will take longer to draw.
#' @param xlim Limits of the horizontal axis.
#' @param ylim Limits of the vertical axis.
#' @param colours Colours of the percentages. Defaults to rainbow colours.
#' @param max.contam Maximum contamination level as a fraction of the sample. Defaults to 0.1 (10\%).
#' @param contam.F14C 14C activity of the sample. Defaults to 'modern' 14C, F14C=1.
#' @param contam.legend Percentages for which numbers will be plotted.
#' @param legend.pos horizontal position beyond which the percentage values will be plotted
#' @param legend.cex font size of the legend
#' @param grid Whether to plot a grid. Defaults to TRUE
#' @param xaxs Whether or not to extend the limits of the horizontal axis. Defaults to \code{xaxs="i"} which does not extend.
#' @param yaxs Whether or not to extend the limits of the vertical axis. Defaults to \code{yaxs="i"} which does not extend.
#' @examples
#'   draw.contamination()
#'   draw.contamination(40e3, 50e3, ka=FALSE)
#' @export
draw.contamination <- function(from=0, to=50e3, ka=TRUE, age.res=500, xlim=c(), ylim=c(), colours=rainbow(age.res), max.contam=.1, contam.F14C=1, contam.legend=max.contam*c(1/100, (1:5)/50, (1:4)/5, 1), legend.pos=.07, legend.cex=0.6, grid=TRUE, xaxs="i", yaxs="i") {
  real.14C <- seq(from, to, length=age.res)
  observed.14C <- seq(0, to, by=diff(real.14C)[1])

  fraction.contaminated <- function(real, observed, F14C=contam.F14C, decimals=5) {
    real.F <- C14toF14C(real, rep(0,length(real)), decimals=decimals)[,1]
    observed.F <- C14toF14C(observed, rep(0,length(real)), decimals=decimals)[,1]
    return((observed.F - real.F ) / (F14C - real.F))
  }
  fractions <- outer(real.14C, observed.14C, fraction.contaminated)

  for(i in 1:length(real.14C))
    for(j in 1:length(observed.14C))
      if(real.14C[i] <= observed.14C[j])
        fractions[i,j] <- NA # only >0 contamination and with F14C>0
  fractions[fractions > max.contam] <- NA # calculate up to a fraction of contamination

  if(ka) {
    ka <- 1e3
    xlab <- expression("real "^14*C~kBP)
    ylab <- expression("observed "^14*C~kBP)
    } else {
      ka <- 1
      xlab <- expression("real "^14*C~BP)
      ylab <- expression("observed "^14*C~BP)
    }

  if(length(xlim) == 0)
    xlim <- extendrange(c(from, to), f=c(0, .1))/ka
  if(length(ylim) == 0)
    ylim <- extendrange(c(to, contaminate(real.14C, 0, max(contam.legend), 0, contam.F14C, 0, MC=FALSE)),
    f=c(0,0.05))/ka

  plot(0, type="n", xlim=xlim, xlab=xlab, ylim=ylim, ylab=ylab, xaxs=xaxs, yaxs=yaxs)
  if(grid)
    grid(lty=2, col=rgb(0,0,0,.2))
  image(real.14C/ka, observed.14C/ka, fractions, col=colours, add=TRUE)
  abline(0, 1, lty=2)

  for(i in contam.legend) {
    contam <- contaminate(real.14C, 0, 100*i, 0, contam.F14C, MC=FALSE)[,1]
    lines(real.14C/ka, contam/ka, lty=3, col=grey(.7))
    text((to+(legend.pos*(to-from)))/ka, max(contam)/ka, labels=paste0(100*i, "%"), cex=legend.cex, adj=c(1,.5))
  }
}
