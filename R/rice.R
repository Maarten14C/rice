# could draw.dates and calibratable be modified to work with custom-built as well as prebuilt curves? E.g., using cc=4 for a custom-built curve, cc=5 for another one...

# show real AMS counts, standards, real samples, multiple measurements

# can intcal.data be made to show the SH data? Currently shows the NH data

# add decay correction delta (difference between year of measurement and year of collection) D = 1000*[(1+D14C/1000)exp(lambda dt) - 1]

# simulate spallation/magnetism, Cs sputtering, AMS beam, stripping

# sample weight functions (per Philippa Ascough's suggestion). Given a %C (perhaps provide estimates for sample types such as peat, bone, ...), a loss during pretreatment, and a required graphite weight, what sample weight will be required?)

# prepare a function to redo deltaR calcs when new Marine curves come out. Using BCADtocalBP(shells$collected), calBPto14C(cc=2) and shells$C14, shells$er. Unclear how the dR errors are obtained.

# fruits-type model that mixes atmospheric and marine calibration curves. Freshwater effects can cause C14 shifts of up to 1k.

# error multipliers, rounding. Could add procedures for different labs, e.g. QUB_bg, etc. This would be useful for reasons of transparency and community standards. Add data from historical UBA standards/backgrounds?



#' @name atom
#' @title Build an atom
#' From the number of protons and neutrons, report the corresponding element/isotope, its abundance, half-life, decay mode and daughter product.
#' @details The number of protons determines the element, and the number of neutrons the isotope. Some combinations of protons and neutrons are impossible (NA), others exist for only a very short time before decaying, and others are stable.
#'
#' If only the number of protons is provided, a range of possible amounts of neutrons is returned, et vice versa.
#'
#' Data downloaded from https://www-nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all
#' @return Invisibly returns a list with the atom's name, symbol, abundance, half-life and decay mode.
#' A message describing the isotope is printed.
#' @param protons The number of protons in the atom's nucleus. Alternatively, the name or symbol of an element can provided.
#' @param neutrons The number of neutrons in the atom's nucleus. If left empty, the function will report the amounts of neutrons of the isotopes of the element in question.
#' @param roundby Rounding of the halflife (if the isotope is unstable). Defaults to 1 decimal.
#' @param talk Report the results. Defaults to TRUE.
#' @author Maarten Blaauw
#' @examples
#'   atom(1,0)
#'   atom("carbon")
#'   atom("C")
#'   atom(6,7)
#'   atom(6,8)
#'   atom("lead")
#'   atom(82, 128)
#' @export
atom <- function(protons=c(), neutrons=c(), roundby=1, talk=TRUE) {

  # if incomplete or wrong entries are provided, tell what more is required
  if(length(protons) == 0 && length(neutrons) == 0) {
    message("please provide a value for protons and/or neutrons, e.g. atom(1,0)")
    return(invisible(NULL))
  }

  if(length(protons) == 1 && is.character(protons)) {
    if(nchar(protons) <= 2)
      idx <- which(tolower(elements$symbol) %in% tolower(protons))[1] else
        idx <- which(tolower(elements$name) %in% tolower(protons))[1]
    if(length(idx) == 0 || is.na(idx)) {
      message("cannot find this element")
      return(invisible(NULL))
    }
    protons <- as.numeric(elements$protons[idx])
  }

  if(length(protons) == 0 && length(neutrons) == 1) {
    if(!neutrons %in% 0:178) {
      message("'neutrons' has to be a whole number between 0 and 178")
      return(invisible(NULL))
    }
    these.protons <- elements$protons[which(elements$neutrons == neutrons)]
    message("for atoms with ", neutrons, " neutrons, choose between ", min(these.protons), 
      " and ", max(these.protons), " protons")
    return(invisible(NULL))
  }

  if(length(protons) == 1 && length(neutrons) == 0) {
    if(!protons %in% 1:118) {
      message("protons should be a whole number between 1 and 118")
      return(invisible(NULL))
    }
    these.neutrons <- elements$neutrons[which(elements$protons == protons)]   
    element.name <- tolower(unique(elements$name[elements$protons==protons]))
    if(length(these.neutrons) == 1)
      message("to see the isotopes for ", element.name, " (", protons, 
        " protons), choose ", these.neutrons, " neutrons") else
          message("isotopes of ", element.name, " (", protons, 
        " protons) have between ", min(these.neutrons), 
        " and ", max(these.neutrons), " neutrons")
    return(invisible(NULL))
  }

  if(length(protons) != 1 || !protons %in% 1:118) {
    message("'protons' has to be a single whole number between 1 and 118")
    return(invisible(NULL))
  }
  if(length(neutrons) != 1 || !neutrons %in% 0:178) {
    message("'neutrons' has to be a whole number between 0 and 178")
    return(invisible(NULL))
  }
    
  # all is OK
  thisone <- which(elements$protons==protons & elements$neutrons==neutrons)
  if(length(thisone) == 0)
    return(NA)

  # convert seconds to human-readable time units
  convert_to_units <- function(s, digits=roundby, seconds_in_day=86400,
    seconds_in_year=365.25*seconds_in_day, seconds_in_millennium=1000*seconds_in_year,
    seconds_in_Ma=1e6*seconds_in_year, seconds_in_Ga=1e9*seconds_in_year) {
    if (s < 1e-9) {
      return(paste0(formatC(s/1e-9, format="e", digits=digits), " nanoseconds"))
    } else if (s < 1e-6) {
      return(paste0(formatC(s/1e-9, format="f", digits=digits), " nanoseconds"))
    } else if (s < 1e-3) {
      return(paste0(formatC(s/1e-6, format="f", digits=digits), " microseconds"))
    } else if (s < 1) {
      return(paste0(formatC(s/1e-3, format="f", digits=digits), " milliseconds"))
    } else if (s < 60) {
      return(paste0(formatC(s, format="f", digits=digits), " seconds"))
    } else if (s < 3600) {
      return(paste0(formatC(s / 60, format="f", digits=digits), " minutes"))
    } else if (s < 100*seconds_in_day) {
      return(paste0(formatC(s / 3600, format="f", digits=digits), " hours"))
    } else if (s < 100*seconds_in_year) {
      return(paste0(formatC(s / seconds_in_day, format="f", digits=digits), " days"))
    } else if (s < 100*seconds_in_millennium) {
      return(paste0(formatC(s / seconds_in_year, format="f", digits=digits), " years"))
    } else if (s < 100*seconds_in_Ma) {
      return(paste0(formatC(s / seconds_in_millennium, format="f", digits=digits), " ka"))
    } else if(s < 100*seconds_in_Ga) {
      return(paste0(formatC(s / seconds_in_Ma, format="f", digits=digits), " Ma"))
    } else if(s / seconds_in_Ga <= 1e4)
      return(paste0(formatC(s / seconds_in_Ga, format="f", digits=digits), " Ga")) else
        return(paste0(formatC(s / seconds_in_Ga, format="e", digits=digits), " Ga"))
  }

  name <- tolower(elements$name[thisone])
  symbol <- elements$symbol[thisone]
  decay <- elements$decay[thisone]
  halflife <- elements$halflife[thisone]
  abundance <- elements$abundance[thisone]
  if(is.na(halflife))
    halflife.txt <- "unknown" else
      if(halflife == 0) # stable
        halflife.txt <- "stable" else
          halflife.txt <- convert_to_units(halflife)

  decay.txt <- "mode not listed"; daughter.txt <- " decay mode not listed"
  new.protons <- c(); new.neutrons <- c()
  decays <- c("stable"="stable", "A"="alpha decay",
    "B+"="beta+ decay", "B-"="beta- decay", 
    "EC"="electron capture",
    "N"="neutron emission", "2N"="two-neutron emission",
    "P"="proton emission", "2P"="two-proton emission",
    "IT"="isomeric transition", "SF"="spontaneous fission",
    "B-N"="beta- and neutron emission",
    "B-2N"="beta- and two-neutron emission",
    "B+P"="beta+ and proton emission",
    "B-A"="beta- and alpha decay",
    "ECP"="electron capture and proton emission",
    "EC+B+"="electron capture / beta-plus decay",
    "2EC"="double electron capture",
    "2B-"="double beta- decay", "2B+"="double beta+ decay",
    "ECSF"="electron capture with spontaneous fission")
  decay.map <- list("stable"=c(0, 0),
    "A"=c(2, 2), "B+"=c(1, -1), "EC"=c(1, -1), "EC+B+"=c(1, -1),
    "B-"=c(-1, 1), "N"=c(0, 1), "2N"=c(0, 2), "P"= c(1, 0),
    "2P"=c(2, 0), "B-N"=c(-1, 2), "B-2N"=c(-1, 3), "B+P"=c(2, -1),
    "ECP"=c(2, -1), "B-A"=c(-1, -3), "2EC"=c(2, -2),
    "2B+"=c(2, -2), "2B-"=c(-2, 2))
  if(!is.na(decay) && decay %in% names(decays))
    decay.txt <- decays[[decay]]

  daughter.protons <- protons
  daughter.neutrons <- neutrons
  daughter.weight <- protons + neutrons
  daughter.name <- NA
  daughter.symbol <- NA
  daughter.txt <- ""
  new.protons <- NA
  new.neutrons <- NA
  if(!is.na(decay)) {
    decay.protons <- 0
    decay.neutrons <- 0
    fission <- NA

    if(decay %in% names(decay.map)) {
      daughter.protons <- daughter.protons - decay.map[[decay]][1]
      daughter.neutrons <- daughter.neutrons - decay.map[[decay]][2]
      daughter.weight <- daughter.protons + daughter.neutrons
    }
    if(decay %in% c("SF", "ECSF")) {
      fission <- " into multiple and variable fragments"
      daughter.protons <- daughter.neutrons <- NA
    }

    if(!decay %in% names(decay.map) && !decay %in% c("SF","ECSF","IT"))
      stop(paste("Unknown decay type:", decay))
  
    if(!is.na(fission))
      daughter.txt <- fission else {
        daughter <- which(elements$protons==daughter.protons & elements$neutrons==daughter.neutrons)
        if(length(daughter) > 0) {
          daughter.name <- tolower(elements$name[daughter])
          daughter.symbol <- elements$symbol[daughter]
          daughter.txt <- paste0(" to ", daughter.name, "-", daughter.weight)
        } else {
            daughter.name <- unique(
              tolower(elements$name[elements$protons == daughter.protons]))
            daughter.symbol <- unique(
              elements$symbol[elements$protons == daughter.protons])
            daughter.txt <- paste0(" to ", daughter.name, "-", daughter.weight)
          }
      }
  }

  if(is.na(halflife))
    msg <- (paste0(name, " (", symbol, "-", protons+neutrons, ", ", 
      protons, " protons and ", neutrons, " neutrons), half-life unknown, ", 
      decay.txt, " to ",
      daughter.name, "-", daughter.weight)) else
        if(halflife == 0)
          msg <- (paste0(name, " (", symbol, "-", protons+neutrons, ", ", 
            protons, " protons and ", neutrons, " neutrons), stable")) else
              msg <- (paste0(name, " (", symbol, "-", protons+neutrons, ", ", 
                protons, " protons and ", neutrons, " neutrons), half-life ",
                  halflife.txt, ", ", decay.txt, daughter.txt))
  if(talk)
    message(msg)

  invisible(list(name=name, symbol=symbol, protons=protons, neutrons=neutrons,
    abundance=abundance, decay=decay,
    halflife=halflife.txt, halflife_seconds=halflife,
    daughter.name=daughter.name, daughter.symbol=daughter.symbol,
    daughter.protons=daughter.protons, daughter.neutrons=daughter.neutrons,
    message=msg))
}



#' @name decay_series
#' @title Decay series of isotopes
#' @description Given a starting isotope, show the decay chain all the way to reaching a stable isotope.
#' @details The `atom` function is called for each of the isotopes in the decay series, and a plot is made. 
#' 
#' @param protons The number of protons in the atom's nucleus. Alternatively, the name or symbol of an element can provided.
#' @param neutrons The number of neutrons in the atom's nucleus. If left empty, the function will report the amounts of neutrons of the isotopes of the element in question.
#' @param plot Plot the results. Defaults to TRUE.
#' @param add Add data to an existing plot (instead of making a new plot). Defaults to FALSE.
#' @param n Number of isotopes to draw, starting from the first isotope. Defaults to drawing all isotopes in the decay series.
#' @param box.size Size of the box around the isotope symbols. Defaults to 0.9. 
#' @param arrows Whether or not to draw arrows between the parent and daughter isotopes. Defaults to TRUE.
#' @param draw.box Whether or not to draw a box around the element. Defaults to TRUE.
#' @param draw.symbol Whether or not to draw the element's symbol, e.g. 14C. Defaults to TRUE.
#' @param grid.col The colour of the grid. Defaults to light grey, \code{grid.col=grey(.75)}.
#' @param xlim Limits of the x axis. Calculated automatically by default.
#' @param ylim Limits of the y axis. Calculated automatically by default.
#' @param cex Controls how much of each isotope box is filled by the isotope's symbol (with superscript mass numbers). Values around 2–3 typically produce labels that occupy most of the box while remaining legible. Defaults to 2.
#' @param offset Offset of the symbol within the isotope's box. Defaults to centre, 0.5.
#' @param mar Margins around the plot. Defaults to mar=c(3,3,0.5,0.5)
#' @param mgp Locations of the axis titles, labels and lines. Defaults to mgp=(1.5, 0.7, 0)
#' @param asp If set to 1, plots the units of the horizontal and vertical axis at 1:1 proportion. 
#' @param talk Report the results. Defaults to TRUE.
#' @author Maarten Blaauw
#' @examples
#'   decay_series("U", 146) # U-238 series
#'   decay_series("U", 143) # U-235 series
#'   decay_series("Th", 142) # Thorium series
#' @export
decay_series <- function(protons, neutrons, plot=TRUE, add=FALSE, n=c(), box.size=.8, arrows=TRUE, draw.box=TRUE, draw.symbol=TRUE, grid.col=grey(.75), xlim=c(), ylim=c(), cex=2, offset=.5, mar=c(3,3,.5, .5), mgp=c(1.5, .7, 0), asp=0, talk=TRUE) {
  series <- list(protons=numeric(), neutrons=numeric(), name=c(), symbol=c(), halflife_s=c(), decay=c(), daughter.protons=c(), daughter.neutrons=c(), message=c())

  snap <- atom(protons, neutrons, talk=FALSE)
  quickplot <- function(p0=snap$protons, n0=snap$neutrons) { # plot the element even if no decay
    if(length(xlim) == 0)
      xlim <- p0*c(.95, 1.05) + c(-1,1) # space for the boxes
    if(length(ylim) == 0)
      ylim <- n0*c(.95, 1.05) + c(-1,1)
    par(bty="l", mar=mar, mgp=mgp)
    if(!add)
      plot(0, type="l", xlim=xlim, ylim=ylim, xlab="protons", ylab="neutrons", asp=asp)
    grid(col=grid.col)
    symbol <- bquote(""^{.(p0+n0)}*.(snap$symbol))
    l.dim <- max(strwidth(symbol), strheight(symbol), box.size) # dimensions
    bx <- box.size/2
    if(draw.box)
      rect(p0, n0, p0+box.size, n0+box.size, col="white", lwd=2, border="black")
    if(draw.symbol)
      text(p0+bx, n0+bx, symbol, cex=cex/(l.dim), adj=offset)
    }

  if(is.null(snap) || length(snap) == 1 && is.na(snap)) {
    if(talk)
      message("no information for this starting element")
    if(plot) quickplot()
    return(invisible(NA))
  }

  if(snap$halflife == "stable") {
    if(talk)
      message("the starting element is stable")
    if(plot) quickplot()
    return(invisible(NA))
  }

  if(!is.na(snap$daughter.protons))
    if(snap$daughter.protons == protons && snap$daughter.neutrons == neutrons) {
      if(talk)
        message("Decay does not change nuclide; terminating chain")
      if(plot) quickplot()
      return(invisible(NA))
    }

  i <- 0
  repeat {
    i <- i+1
    series$name[[i]] <- snap$name
    series$symbol[[i]] <- snap$symbol
    if(snap$halflife == "stable")
      snap$halflife <- NA
    series$halflife[[i]] <- snap$halflife
    series$halflife_s[[i]] <- snap$halflife_seconds
    series$message[[i]] <- snap$message
    series$protons[[i]] <- snap$protons
    series$neutrons[[i]] <- snap$neutrons
    next.p <- snap$daughter.protons
    next.n <- snap$daughter.neutrons
    series$daughter.protons[[i]] <- next.p
    series$daughter.neutrons[[i]] <- next.n
    series$decay[[i]] <- snap$decay

    if(is.na(snap$daughter.protons) || is.na(snap$daughter.neutrons))
      break
    if(is.na(snap$decay) || snap$decay == "stable")
      break
    if(snap$decay == "IT") # the daughter isotope will be the same as the parent
      break
    if(snap$daughter.protons == snap$protons &&
      snap$daughter.neutrons == snap$neutrons) {
        if(talk)
          message("Decay does not lead to a unique daughter nuclide; terminating chain")
        break
    }

    snap <- atom(next.p, next.n, talk=FALSE)

    if(length(snap) == 1 && is.na(snap)) {
      series$name[[i+1]] <- paste0(
        unique(elements$name[elements$protons == next.p]))
      series$symbol[[i+1]] <- unique(
        elements$symbol[elements$protons == next.p])
      series$protons[[i+1]] <- next.p
      series$neutrons[[i+1]] <- next.n
      series$halflife[[i+1]] <- NA
      series$message[[i+1]] <- paste0(tolower(series$name[[i+1]]),
        "-", next.p+next.n, " (not in database)")
      break
    }
  }

  if(talk) {
    msg <- c()
    for(i in seq_along(series$protons)) {
      decay.txt <- strsplit(series$message[[i]], " to ")[[1]][1]
      msg[i] <- paste0(strrep(" ", i), decay.txt,
        if(i<length(series$message)){" to:\n"}) # add an extra space for each step
      }
    message(msg)
    }
  
  if(plot) {
    rng.protons <- c(.95, 1.05)*
      range(c(series$protons, series$daughter.protons), na.rm=TRUE)
    rng.neutrons <- c(.95, 1.05)*
      range(c(series$neutrons, series$daughter.neutrons), na.rm=TRUE)
  
    decay.cols <- c("B-"="blue", "EC+B+"="orange", A="red", 
     EC="grey", P="darkgreen", SF="brown", N="cyan")

    bx <- box.size/2
    if(length(xlim) == 0)
      xlim <- range(rng.protons+c(-1,1)) # space for the boxes
    if(length(ylim) == 0)
      ylim <- range(rng.neutrons+c(-1,1))
    par(bty="l", mar=mar, mgp=mgp)
    if(!add)
      plot(0, type="l", xlim=xlim, ylim=ylim, xlab="protons", ylab="neutrons", asp=asp)  
   grid(col=grid.col)
   if(length(n)==0)
      n <- length(series$protons)
    for(i in 1:min(n, length(series$decay))) {
      p <- series$protons[[i]]
      n0 <- series$neutrons[[i]]
      hl <- series$halflife[[i]]
      dp <- series$daughter.protons[[i]]
      dn <- series$daughter.neutrons[[i]]

      d <- if(i <= length(series$decay)) series$decay[[i]] else NA
      col <- 1
      if(!is.na(hl) && hl != "stable")
        if(!is.na(d) && d %in% names(decay.cols))
         col <- decay.cols[[d]]
      if(arrows &&
        !anyNA(c(p, n0, dp, dn)) && (p != dp || n0 != dn))
          arrows(p+bx, n0+bx, dp+bx, dn+bx, col=col, code=2, length=.05)
    }

    nseq <- seq_len(min(n, length(series$protons)))
    l.dim <- 0
    for(i in nseq) {
      np <- series$protons[[i]] + series$neutrons[[i]]
      symbol <- bquote(""^{.(np)}*series$symbol[[i]]) #, e.g., ^{14}C
      l.dim <- max(l.dim, strwidth(symbol), strheight(symbol)) # dimensions
    }

    for(i in nseq) {
      p <- series$protons[[i]]
      n0 <- series$neutrons[[i]]
      hl <- series$halflife[[i]]
      col <- 1 
      legend.cols <- col
      if(hl != "stable" && !is.na(hl)) {
        d <- if(i <= length(series$decay)) series$decay[[i]] else NA
        if(d %in% names(decay.cols))
          col <- decay.cols[[d]]
          legend.cols[i] <- col   
      }
      if(draw.box)
        rect(p, n0, p+box.size, n0+box.size, col="white", lwd=ifelse(i==1,2,1), border=col)
      if(draw.symbol)
        text(p+bx, n0+bx, bquote(""^{.(p+n0)}*.(series$symbol[[i]])),
          cex=cex/(l.dim*box.size), adj=offset)
    }

    decay.legend <- unique(unlist(series$decay))
    decay.labels <- sapply(decay.legend, function(x) # alpha/beta symbols
      switch(x, "A"=expression(alpha), "B-"=expression(beta^"-"), 
        "B+"=expression(beta^"+"), "EC+B+"=expression("EC+"*beta^"+"),
        "B-N"=expression(beta^"-"*N),"B-2N"=expression(beta^"-"*"2N"),
        "B+P"=expression(beta^"+"*"P"), "B-A"=expression(beta^"-"*alpha),
         "2B+"=expression(2*beta^"+"), "2B-"=expression(2*beta^"-"),  x))
    legend.cols <- unname(decay.cols[decay.legend])
    legend("topleft", legend=decay.labels, text.col=legend.cols, bty="n", bg="white") 
  }
  
  invisible(series)
}



#' @name howmuchC14
#' @title Amount of C14 particles in a sample
#' @description Calculate the expected amount of remaining C14 atoms in a sample, given its weight and age.
#' @details The number of carbon atoms in the sample is estimated. Given the known C14/C ratio at F=1, and given the sample's age, we can estimate the number of remaining C14 atoms. Given a 12C current at the detector end of an AMS, we can then also calculate how many 14C ions would be counted per second and minute. 
#' Measured C14 activities/concentrations are expressed relative to 95% of the 14C/C ratio of the Oxalic Acid I (OxI) standard (and normalised to d13C=-25 permille). This standard was made of oxalic acid from a 1957 beet sugar crop, distributed by the US National Bureau of Standards. When this standard ran out, in 1977 a French beetroot crop was prepared the same way and distributed as OxII (cross-calibrated to OxI and still in use by most C-14 laboratories).  
#' This standard ratio, "modern" (F14C=1) is 1.176e-12, or in other words, c. 1.2 C-14 atoms per 1 trillion carbon atoms. Note that this is a conventional reference value, not the actual atmospheric 14C/C ratio in AD 1950. Due to the Suess effect (dilution by CO2 from C14-free fossil fuel combustion), the real atmospheric ratio in AD 1950 was already slightly lower.
#'
#' Note that backgrounds are not modelled (but could be investigated by e.g. typing \code{howmuchC14(45e3)} which gives on the order of 1 background count per second). The calculated C14 count rate assumes no isotopic fractionation.
#' @return The estimated number of C14 atoms.
#' @param age The age of the sample (in cal BP per default, or in C14 BP if use.cc=FALSE).
#' @param is.F By default, ages are assumed to be in either cal BP or 14C BP. If \code{is.F=TRUE}, age is assumed to be on the F14C scale.
#' @param wght The weight of the sample (in mg). Defaults to 1 mg.
#' @param use.cc Whether or not to use the calibration curve. If set to \code{use.cc=FALSE}, then we assume that the age is the radiocarbon age (this enables ages beyond the reach of the calibration curves to be used).
#' @param Av Avogadro's number, used to calculate the number of carbon atoms in the sample.
#' @param C14.1950 The standard 14C/C ratio at 0 cal BP (AD 1950), defined as 95\% of the C-14 activity of NBS Oxalic Acid I (oxI), normalized to d13C=–25 permille. 
#' @param current The beam current of 12C+ ions as measured at the Faraday cup (C12 detector). Defaults to \code{current=25e-6}, 25 microamperes, a typical value for graphite targets on modern AMS systems. Gas targets generally yield c. 4-5 times lower currents.
#' @param format The format of the printed numbers. Defaults to either scientific (for large numbers) or as fixed-point, depending on the size of the number.
#' @param cc calibration curve for C14 (see \code{caldist()}).
#' @param postbomb Whether or not to use a postbomb curve (see \code{caldist()}).
#' @param glue Glue postbomb and prebomb curves together. Defaults to 0 (none), can be 1 (IntCal20 + NH1), 2 (IntCal20 + NH2), 3 (IntCal20 + NH3), 4 (SHCal20 + SH1-2) or 5 (SHCal20 + SH3). Note that this will override the value of cc. Can also be used to glue other pre-bomb and post-bomb curves, e.g., \code{glue="NH1_monthly"}.
#' @param cc.dir Directory of the calibration curves. Defaults to where the package's files are stored (system.file), but can be set to, e.g., \code{cc.dir="curves"}.
#' @param thiscurve As an alternative to providing cc and/or postbomb, the data of a specific curve can be provided (3 columns: cal BP, C14 age, error).
#' @param talk Whether or not to provide feedback (defaults to TRUE).
#' @param as.AMS If set to true, will calculate how many atoms would be counted in an AMS 
#' @param decimals Number of decimals to be returned for F and atom counts.
#' @author Maarten Blaauw
#' @examples
#'   howmuchC14(0) # recent sample
#'   howmuchC14(55e3) # at dating limit
#'   howmuchC14(145e3) # way beyond the dating limit, 1 C14 atom per mg remains
#' @export
howmuchC14 <- function(age, wght=1, is.F=FALSE, use.cc=TRUE, Av=6.02214076e23, C14.1950=1.176e-12, current=25e-6, format="g", cc=1, postbomb=FALSE, glue=0, cc.dir=NULL, thiscurve=NULL, talk=TRUE, as.AMS=TRUE, decimals=3) {
  if(length(age) > 1)
    stop("can only handle single value for age")
  
  if(is.F)
    F <- age else {
      if(use.cc) {
         F <- calBPtoF14C(age, cc=cc, postbomb=postbomb, glue=glue, cc.dir=cc.dir, thiscurve=thiscurve)[,1]
         if(is.na(F)) {
           message("Cannot use calibration curve for this age, assuming C14 age")
           F <- C14toF14C(age)
         }
       } else
           F <- C14toF14C(age) # then 'age' is on the C14 scale
    }
  F <- as.numeric(F)
  
  atoms <- (wght/1e3)*Av/12 # number of C atoms in a mg
  C14 <- as.numeric(ceiling(F * C14.1950 * atoms)) # C14 atoms rounded up
  
  # from the current, calculate amount of ions ending up at the C12 Faraday cup:
  i12 <- current/1.602e-19 # e, electric charge for a single proton, in coulomb
  # from this, we can model the numbers of C14 atoms arriving at the C14 detector:
  i14 <- C14.1950 * i12 * F # charge in A = protons arriving per second
  persecond <- round(i14, decimals)

  atoms <- formatC(atoms, format=format, digits=decimals)
  C14.talk <- formatC(C14, format=format, digits=decimals)
  decays_val <- C14 * log(2) / (5730 * 365.25) # per day
  decays_str <- formatC(decays_val, format="fg", digits=decimals)

  if(talk) {
    message(wght, " mg carbon contains ", atoms, " C atoms")
    if(is.F) 
      message("14C atoms remaining at ", age, " F14C: ", C14.talk) else {
        if(use.cc)
          message("14C atoms remaining at ", round(age, decimals), " cal BP (F=", round(F, decimals), "): ", C14.talk) else
            message("14C atoms remaining at ", round(age, decimals), " 14C BP (F=", round(F, decimals), "): ", C14.talk)
        }
    message(decays_str, " 14C atoms in the sample will decay each day")
    if(as.AMS)
      message(paste0("For a 12C current of ", current*1e6, " micro-ampere (", formatC(i12, format=format, digits=decimals),
        " 12C/second) at the AMS detector,\n  ",
        persecond, " 14C particles would be counted per second (",
        60*persecond, " per minute)" ))
  }

  invisible(list(C14=C14, C14persecond=i14, decaysperday=decays_val, i12=i12, C14.1950=C14.1950, totC=as.numeric(atoms)))
}



#' @name radio
#' @title Listen to radiocarbon being measured
#' @description A sonification of the amount of C14 being detected in an AMS. 
#' @details The sonification can be heard as 'Geiger counter-like' clicks (individual C14 atoms ending up in the detector) and/or as a possibly slightly meandering sine wave (e.g., for a sample of age 5000 14C BP, we'd count c. 105.6 atoms per second in the AMS, so this would become a tone of frequency 105.6 Hz). The calculations are based on the function \code{howmuchC14}.
#'
#' Instead of counting the particles themselves (as done in an AMS), we can also count the decay events as used in radiometric counting. In this case, we need much larger samples (wght) and longer counting times (duration).
#'
#' To produce and listen to the sounds, the package 'audio' (and possibly 'tuneR') has to be installed. A warning will be provided if it isn't installed.
#' Note that the calculations in this function do not include error multipliers or corrections such as for fractionation and backgrounds. 
#' @return A sound (tone and or clicks) is played and returned invisibly.
#' @param age The age of the sample (in C14 BP by default, or in cal BP if use.cc=FALSE).
#' @param duration How long the sample will count (and sound) for in seconds. Defaults to 10 seconds.
#' @param duration.unit Unit of the duration. Left empty by default (c()) which then makes seconds the unit ('s'), or if as.decays is used, days ('d'). Even if set to 'd', the sound will be played back using the seconds unit ('s'). 
#' @param is.F By default, ages are assumed to be in either cal BP or 14C BP. If \code{is.F=TRUE}, age is assumed to be on the F14C scale.
#' @param use.cc Whether or not to use the calibration curve. If set to \code{use.cc=FALSE} (the default), then we assume that the age is the radiocarbon age (this enables ages beyond the reach of the calibration curves to be used).
#' @param as.decays Work with the C14 decays. Defaults to FALSE, which works with the number of C14 atoms instead. If set to true, you'll need larger sample sizes (wght) and longer counting times (duration) to get decent counts.
#' @param wght The weight of the sample (in mg). Defaults to 1 mg.
#' @param play Whether or not to play the sound
#' @param as.clicks Make the C-14 counts sound as clicks, based on random Poisson sampling. Defaults to TRUE.
#' @param click_length Length of the clicks. Defaults to 80, but can be altered to make the clicks sound different (larger values leave a larger 'echo')
#' @param as.tone Make the C-14 frequency (counts per second) sound as a wave (e.g., 105 cps becomes a 105 Hz sine wave). This can either be a constant wave or be meandering (see `wobble`).
#' @param tone.volume Volume of the tone/wave relative to that of the clicks.
#' @param wobble Drift of the tone along the mean. Defaults to 0.00001. Increasing this value can cause the sound to change from a stable tone (low values of `wobble`) via an unstable one to noise and eventually clicks. Values outside the 0-1 range are unlikely to work as expected. 
#' @param sr Sampling rate. This audio quality defaults to 44100 (per second), which is based on the CD standard.
#' @param visualise the counts and calculation of the C14 age. Defaults to TRUE.
#' @param cex Size of the font. Defaults to \code{cex=0.6}.
#' @param return.sound Return the sound as an invisible object. If set to FALSE (the default), instead returns the counts and the calculated C14 age.
#' @param ... Optional constants to be entered into the function `howmuchC14`
#' @author Maarten Blaauw
#' @examples
#'   radio(0)
#'   radio(45000)
#'   # decay events over 1 hour in 1 gram of carbon of age 500 14C BP:
#'   radio(500, as.decays=TRUE, wght=1000, duration=1/24) 
#' @export
radio <- function(age, duration=10, duration.unit=c(), is.F=FALSE, use.cc=FALSE, as.decays=FALSE, wght=1, play=interactive(), as.clicks=TRUE, click_length=80, as.tone=TRUE, tone.volume=0.5, wobble=c(), sr=44100, visualise=TRUE, cex=.5, return.sound=FALSE, ...) {

  hasaudio <- requireNamespace("audio", quietly=TRUE)
  hastuneR <- requireNamespace("tuneR", quietly=TRUE)
  if(!hasaudio && !hastuneR) {
    if(interactive())
      install.packages(c("audio", "tuneR")) else
        stop("Please install the audio or tuneR packages:\ninstall.packages(c(\"audio\", \"tuneR\"))")
  }

  if(as.decays && as.tone && play)
    message("Tone not supported in as.decay mode")

  if(!as.clicks && !as.tone && play)
    message("No sound to play")

  if(length(duration.unit) == 0)
    duration.unit <- if(as.decays) "d" else "s"
  if(as.decays)
    as.tone <- FALSE

  counts <- howmuchC14(age, wght=wght, as.AMS=!as.decays, use.cc=use.cc, is.F=is.F, ...)

  if(as.decays) # rate = per second
    rate <- counts$decaysperday / (60*60*24) else # C14 decays
      rate <- counts$C14persecond # C14 atoms
  duration_sec <- duration * if(duration.unit == "d") 86400 else 1

  C14_count <- rpois(1, rate * duration_sec) # simulate amount of particles...
  event_times <- sort(runif(C14_count, 0, duration_sec)) # and their timing
  n <- sr * duration_sec # number of steps

  if(visualise) {
    if(as.decays) { # radiometric: counting decays, not atoms
      A <- C14_count / duration
      if(duration.unit == "d") # decays per day
        A0 <- counts$decaysperday else
          A0 <- counts$decaysperday / (24*3600) # per second
      asF <- A / A0

      has_counts <- C14_count > 0
      if(has_counts)
        asF.sd <- asF / sqrt(C14_count) else
          asF.sd <- NA
      as.C14 <- F14CtoC14(asF, if(has_counts) asF.sd, roundby=0)
    } else { # AMS: counting atoms
        C12_count <- counts$i12 * duration
        AN <- C14_count / C12_count
        asF <- AN / counts$C14.1950
        has_counts <- C14_count > 0
        if(has_counts) {
          AN.sd <- AN * sqrt(1/C14_count + 1/C12_count)
          asF.sd <- asF / sqrt(C14_count)
      } else {
          AN.sd <- NA
          asF.sd <- NA
        }
        as.C14 <- F14CtoC14(asF, if(has_counts) asF.sd, roundby=0)
      }

   pm <- function(x, sd, sci=FALSE) # print x +- y
     paste(format(x, format="fg", scientific=sci), "\u00B1", format(sd, scientific=sci))

    if(as.decays) {
      labels <- c(expression({}^{14}*C ~ decays ~ (A %+-% sqrt(A))),
        expression(relative~to~AD~1950~(F == A / A[0])),
        expression(as ~ {}^{14}*C ~ BP == -8033 %.% ln(F)))
      } else {
        labels <- c(
          expression({}^{14}*C ~ counts ~ (N %+-% sqrt(N))),
          expression({}^{12}*C ~ counts ~ (N %+-% sqrt(N))),
          expression({}^{14}*C/{}^{12}*C ~ ratio ~ (A[N])),
          expression(relative ~ to ~ AD ~ 1950 ~ (F == A[N]/A[0])),
          expression(as ~ {}^{14}*C ~ BP == -8033 %.% ln(F)))
        }

    values <- if(as.decays)
      c(pm(round(C14_count), round(sqrt(C14_count))), pm(asF, asF.sd), paste(as.C14[1], "\u00B1", as.C14[2])) else
        c(pm(round(C14_count), round(sqrt(C14_count))), pm(round(C12_count), round(sqrt(C12_count)), TRUE),
          pm(AN, AN.sd, TRUE), pm(asF, asF.sd), paste(as.C14[1], "\u00B1", as.C14[2]))

    layout(matrix(c(1,2), nrow=2)) # first the equations
    plot(0, type="n", xlim=c(0,1), ylim=c(0,1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
    nl <- length(labels)
    cex_val <- par("pin")[2] / (2*nl)
    cex_val <- min(1, max(0.87, cex_val))
    line_height <- strheight("M", cex = cex_val) * 1.7
    y <- .95 - (0:(nl-1)) * line_height

    if(as.decays) {
      y <- y[1:3]
      cols <- c(1, grey(0.5), 1)
    } else
        cols <- c(1, rep(grey(0.5), 3), 1)

    text(0.02, y, labels, cex=cex_val, adj=0, col=cols)
    text(0.98, y, values, cex=cex_val, adj=1, col=cols)

    par(mar=c(3, 1, 1, 1), mgp=c(2, .7,0))
    xlab <- ifelse(duration.unit == "d", "duration (d)", "duration (s)")
    plot(0, type="n", xlim=c(0, duration), ylim=c(0, 1), xlab=xlab, ylab="", yaxt="n", bty="n")
    alpha <- 1/(1+((length(event_times))/200)^0.95) # transparency of the lines

    event_times_plot <- if(duration.unit == "d")
      event_times / (3600 * 24) else
        event_times
    if(length(event_times_plot) > 0)
      segments(x0=event_times_plot, y0=0, x1=event_times_plot, y1=1, lwd=.5, col=rgb(0,0,0,alpha))
  }

  # make random clicks based on the rate
  if(as.clicks) {
    audio_duration <- duration  # e.g. 10 seconds of playback
    n <- sr * audio_duration # number of bins
    scale <- audio_duration / duration_sec # if in days, replay in seconds
    event_times_audio <- event_times * scale # times of events, in seconds
    event_idx <- floor(event_times_audio * sr) + 1 # time bin of event

    click <- runif(1, .5, 1.5) * rnorm(click_length) * exp(-seq(0, 6, length.out=click_length)) # decaying noise
    clicks <- numeric(n)
    for(i in event_idx) { # when a click happens...
      end <- min(i + click_length-1, n) # find its time bins...
      clicks[i:end] <- clicks[i:end] + click[1:(end-i+1)] # superpose it to any other clicks at that time
    }
    if(max(abs(clicks), na.rm=TRUE) > 0)
      clicks <- clicks / max(abs(clicks), na.rm=TRUE) # normalise
  }

  # turn the count rate (counts/s=Hz) into a tone (will not be heard at low count rates)
  if(as.tone) {
    if(is.null(wobble)) # default perceptual smoothing
      wobble <- rate / (sr * 2000) # frequency of events relative to sampling rate  
    rate_est <- numeric(n)
    counts <- tabulate(floor(event_times * sr)+1, nbins=n) # find the time bins
    rate_est[1] <- rate # start the tone at frequency of 'rate'
    for(i in 2:n) {
      inst_rate <- counts[i] * sr # find the rate at time i
      rate_est[i] <- (1-wobble)*rate_est[i-1] + wobble*inst_rate # mean-reverting stochastic process
    }
    tone <- sin(cumsum(2 * pi * rate_est / sr)) # wobbly wave
  }
  
  if(as.clicks && as.tone)
    sound <- rbind(tone.volume*tone, clicks) else
      if(as.clicks)
        sound <- clicks else 
          sound <- tone

  if(play && (as.clicks+as.tone>0)) {
   sys <- Sys.info()[["sysname"]] # Windows, Mac or Linux
   if(sys %in% c("Darwin", "Windows") && nrow(audio::audio.drivers()) > 0) # can use 'play'
     audio::play(audio::audioSample(sound, sr)) else {
      if(is.matrix(sound))
        sound <- sound[1,] + sound[2,] # reformat sound
      mx <- max(abs(sound))
      if(mx > 0)
        sound <- sound / mx

      wave <- tuneR::Wave(left=as.integer(sound * 32767), samp.rate=sr, bit=16) # 32767 = 16 bit
      tmpwav <- tempfile(fileext=".wav")
      tuneR::writeWave(wave, tmpwav)
      player <- Sys.which(c("ffplay","pw-play","paplay","aplay"))
      player <- player[nzchar(player)][1]
      if(is.na(player))
        stop("No audio player found")
      system(paste(player, "-nodisp -vn -autoexit -loglevel quiet -hide_banner",
        shQuote(tmpwav), "> /dev/null 2>&1"))
     }
  }

  if(return.sound)
    invisible(sound) else
      invisible(list(C14_count=C14_count, duration=duration, F14C=asF, C14_age=as.C14))
}



#' @name C14.cycle
#' @title Fill the atmosphere and ocean with 14C
#' @description Simulates the accumulation and exchange of 14C within the atmosphere and ocean over time, based on annual production, decay and exchange between the atmosphere and ocean.
#' @details This is a very simplistic 'toy' model with two layers (atmosphere and ocean) between which carbon is exchanged. The model ignores essential carbon cycling features such as vegetation, soils, ocean layers, circulation and upwelling and should not be taken as anything other than an educational tool. 
#'
#' The simulation starts with no C14 present in either of the reservoirs (but both reservoirs will have plenty C12). As time passes, C14 produced by cosmic rays accumulates inside the atmosphere, with some leaking to the ocean at each step. Also at each step, some of the C14 will decay (at a rate dictated by its halflife). At each step, the total mass of 14C in both reservoirs is returned as well as its proportion as 14C/12C. Users can play around with values of annual 14C production, step duration, total mass, 12C concentration, exchange between the atmosphere and ocean, and half-life.
#'
#' Estimates of how much C14 is produced per year vary, but can be approached as follows. Up in the atmosphere (production peaks at c. 10 km height, where airplanes fly), cosmic ray collisions produce c. 1-2 14C atoms per second per cm^2 on average (c. 17k/m2; actual rates vary and are higher near the poles). Since the surface of a sphere equals 4*pi*r^2 and given Earth's radius of 6371 km, the amount of C14 particles produced per year would be n = 17e3 * 365.25 * 24*60*60 * 4*pi*(1e3*(6371+10))^2 or c. 2.745e26 C14 particles. Multiplied by the mass of a single C14 atom, n * 2.32e-26 = c. 6.37 kg. We round this to 6 kg, about the weight of a well-fed (i.e., slightly chonky) cat. 
#' 
#' Default values for the total masses of the atmosphere and the ocean are from Wikipedia, and C-12 concentrations are set at 300 ppm (i.e., pre-industrial atmospheric CO2 concentrations). With these values and the above C-14 production rate, the atmosphere contains around 480 kg of C14 at equilibrium, about the mass of a polar bear. A further 50,000 kg of C-14 resides in the ocean, roughly the mass of 10 African Elephants. The coefficient that sets the rate of exchange between the atmosphere and ocean is tuned such that with default settings, 14C/12C ratios reach expected values of c. 1.2e-12. 
#' @return A plot and the underlying values of the atmospheric and ocean 14C concentration (in kg and as proportion of 12C).
#' @param rate Annual production of 14C, in kg. Defaults to 6 kg/yr (see details). Instead of a constant, can also be provided as a time-series (e.g., a random walk). 
#' @param duration Amount of time over which to cycle, in years.
#' @param n.steps Number of steps.
#' @param halflife Half-life of 14C. Defaults to Libby Halflife of 5730 yr, but can also be set to much shorter or longer to see the impacts on concentrations.
#' @param CO2.ppm.atm Concentration (in ppm) of atmospheric 12C (as CO2). Defaults to 300 ppm. 
#' @param CO2.ppm.ocean Concentration (in ppm) of surface ocean 12C. Defaults to 300 ppm. 
#' @param mass.atm Mass (in kg) of carbon in the atmosphere (as CO2).
#' @param mass.ocean Mass (in kg) of carbon in the ocean.
#' @param exchange Exchange rate of carbon (as CO2) between the atmosphere and the ocean. Defaults to 0.0195 (since this results in reasonable atmospheric 14C concentrations).
#' @param f.ocean Fraction of the ocean partaking in exchange with the atmosphere. Implicitly models surface and deep ocean. Defaults to 0.39.
#' @param as.ratio Plot carbon in the reservoirs as 14C/12C ratio (default). If set to FALSE, is plotted as mass of 14C.
#' @param C14.1950 The standard 14C/C ratio at 0 cal BP (AD 1950), defined as 95\% of the C-14 activity of NBS Oxalic Acid I (oxI), normalized to d13C=–25 permille. Set to NA to avoid plotting the dashed horizontal line. 
#' @param col.atmosphere Colour of the curve depicting the atmosphere's values. Defaults to sunny orange.
#' @param col.ocean Colour of the curve depicting the ocean's values. Defaults to blue.
#' @param col.rate Colour of the rate (annual addition of C14 to the atmosphere). Defaults to dark grey.
#' @param bg.rate Colour of the panel behind the rate curve (annual addition of C14 to the atmosphere). Defaults to light grey.
#' @param x.lim Axis limits of the horizontal axis. Calculated automatically by default.
#' @param y.lim Axis limits of the vertical axis. Calculated automatically by default.
#' @author Maarten Blaauw
#' @examples
#'   C14.cycle()
#'   mn <- 6; rw <- mn; set.seed(67)
#'   for(i in 2:1e5) rw[i] <- rw[i-1] + rnorm(1, 0, .001)
#'   C14.cycle(rw)
#' @export
C14.cycle <- function(rate=6, duration=50e3, n.steps=1e5, halflife=5730, CO2.ppm.atm=300, CO2.ppm.ocean=300, mass.atm=5e18, mass.ocean=1.4e21, exchange=0.0195, f.ocean=0.39, as.ratio=TRUE, C14.1950=1.176e-12, col.atmosphere="orange", col.ocean="blue", col.rate="darkgrey", bg.rate=rgb(0,0,0,.05), x.lim=c(), y.lim=c()) {

  # total mass of C (mostly 12C) in the two reservoirs
  C12.atm <- mass.atm * CO2.ppm.atm * 1e-6 * (12/44) 
  C12.ocean <- mass.ocean * CO2.ppm.ocean * 1e-6 * (12/44)

  if(duration<100e3)
    n.steps <- max(2, duration) # then take yearly steps (at least two)
  if(length(rate) > 1) # production rate can vary over time
    rate <- approx(1:length(rate), rate, 1:n.steps)$y  
  
  # time vector (years)
  time <- seq(0, duration, length=n.steps)
  stepsize <- time[2]-time[1]

  # 14C storage vectors. Initial 14C = 0
  atm <- ocean <- numeric(n.steps)
  atm[1] <- ocean[1] <- 0

  for(i in 2:(n.steps)) {
    decay.atm <- log(2)/halflife * atm[i-1]
    decay.ocean <- log(2)/halflife * ocean[i-1]

    R.atm <- atm[i-1] / C12.atm # previous 14C/12C ratio in the atmosphere...
    R.ocean <- ocean[i-1] / C12.ocean # and in the ocean

    # fluxes move carbon proportional to ratio differences. 'exchange' steers how fast this exchange goes    
    flux.atm.ocean <- exchange * (R.atm - R.ocean) * C12.atm
    prod <- if(length(rate)==1) rate else rate[i]
    atm[i] <- atm[i-1] + (prod - decay.atm - flux.atm.ocean) * stepsize # fill in this step
    ocean[i] <- ocean[i-1] + (flux.atm.ocean - decay.ocean) * stepsize # fill in this step
  }

  R.atm <- atm / C12.atm
  R.ocean <- ocean / (f.ocean * C12.ocean) # not all of the ocean interacts with the atmosphere
  
  if(length(x.lim) == 0)
    xlim <- range(time)   
  
  if(as.ratio) {
    if(length(y.lim) == 0)
      y.lim <- c(0, max(R.atm, R.ocean))  
    plot(time, R.atm, type="l", col=col.atmosphere, xlab="time (yr)", ylab="14C/12C ratio", xlim=x.lim, ylim=y.lim)
    lines(time, R.ocean, col=col.ocean)
    abline(h=C14.1950, col=col.atmosphere, lty=2)
  } else {
      if(length(y.lim) == 0)
        y.lim <- c(0, max(atm, ocean))  
      plot(time, atm, type="l", col=col.atmosphere, xlab="time (yr)", ylab="14C (kg)", xlim=x.lim, ylim=y.lim)
      lines(time, ocean, col=col.ocean)
      abline(h=C14.1950*mass.atm*CO2.ppm.atm/1e6 * (12/44), col=col.atmosphere, lty=2)
    }

  if(length(rate) > 1) {
    yr <- par('usr')
    rr <- range(rate)  
    if(rr[1] == rr[2]) # rate is constant
      x <- rep(mean(yr), length(rate)) else
        x <- yr[3] + (rate - rr[1]) * diff(yr[3:4]) / diff(rr) # scale to lefthand axis
    rect(yr[1], yr[3], yr[2], .35*yr[4], col=bg.rate, border=NA)
    lines(time, .3*x, col=col.rate) # only go up to one third the height of the lefthand axis
    axis(4, at=c(0, .3*(yr[4]-yr[3])), labels=c(round(min(rate),2), round(max(rate),2)), cex.axis=.8, col.axis=col.rate, col.ticks=col.rate, padj=-1.2)
    legend("right", legend=c("atmosphere", "ocean", "rate"),
      col=c(col.atmosphere, col.ocean, col.rate), lty=1, bty="n", cex=0.7)
  } else 
    legend("right", legend=c("atmosphere", "ocean"),
      col=c(col.atmosphere, col.ocean), lty=1, bty="n", cex=0.7)

  invisible(data.frame(time=time, rate=rate, atm=atm, ocean=ocean, R.atm=R.atm, R.ocean=R.ocean))
}



#' @name adjust.fractionation
#' @title Adjust a radiocarbon age for fractionation
#' @description Calculate the radiocarbon age by adjusting a sample's d13C to the reference d13C of -25 permil. It is planned to update this function to more properly reflect calculations in the 14CHRONO lab.
#' @details Radiocarbon ages are corrected for fractionation (which can take place in the field, or during lab pretreatment and measurement), by calculating the radiocarbon age as if the d13C fractionation were at the d13C of the standard (-25 permil). Errors are not taken into account.
#' @return The fractionation-adjusted age.
#' @param y The age of the sample (in C14 by default, but can also be in F or pMC).
#' @param d13C The measured d13C value.
#' @param reference_d13C The reference/standard d13C value (OX2, oxalic acid 2, NIST SRM 4990C made from 1977 French beet molasses), set at -25 permil by default.
#' @param timescale Type of radiocarbon age. Can be in `C14` (default), `F14C` or `pMC`.
#' @author Maarten Blaauw
#' @examples
#'   adjust.fractionation(5000, -17)
#' @export
adjust.fractionation <- function(y, d13C, reference_d13C=-25, timescale="C14") {
  ratio = (1 + reference_d13C / 1000) / (1 + d13C / 1000)
  timescale <- tolower(timescale)

  if(grepl("^c", timescale))
    return(-8033 * log( ratio^2 * exp(-y / 8033) ))
  if(grepl("^f", timescale))
    return(y * ratio^2)
  if(grepl("^p", timescale))
    return((y/100) * ratio^2 * 100) # calculate as F

  stop("Unknown timescale; use 'C14', 'F', or 'pMC'")
}



#' @name adjust.background
#' @title Adjust a radiocarbon age for background measurements
#' @description Calculate the radiocarbon age by adjusting it for a measured background. It is planned to update this function to more properly reflect calculations in the 14CHRONO lab.
#' @details Radiocarbon ages are measured using a series of standards and backgrounds, and the raw values are then corrected for these background values. Backgrounds are >0 (in F14C) owing to contamination in even the cleanest lab.
#' @return The background-adjusted age.
#' @param y The age of the sample (in C14 by default, but can also be in F or pMC).
#' @param er The error of the date.
#' @param bg The background measurement. Should be in the same timescale as that of the sample.
#' @param bg.er The error of the background measurement. Should be in the same timescale as that of the sample.
#' @param timescale Type of radiocarbon age. Can be in `C14` (default), `F14C` or `pMC`.
#' @author Maarten Blaauw
#' @examples
#'   adjust.background(9000, 50, 45000, 200)
#' @export
adjust.background <- function(y, er, bg, bg.er, timescale="C14") {
  timescale <- tolower(timescale)
  is_c14 <- grepl("^c", timescale)
  is_pmc <- grepl("^p", timescale)
  is_f14c <- grepl("^f", timescale)

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
      stop("Unknown timescale; use 'C14', 'F', or 'pMC'")

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



# to read in the required calibration curve(s)
build.curve <- function(y=c(), er=c(), cc=1, thiscurve=c(), cc.dir=c(), is.F=FALSE, is.pMC=FALSE, postbomb=0, glue=0, bombalert=TRUE, cc.resample=FALSE, cc0.res=1e3) {
  if(length(thiscurve) > 0)
    return(thiscurve) else {
      if(cc == 0) { # no ccurve needed
        xseq <- seq(y-4*er, y+4*er, length=cc0.res)
        return(cbind(xseq, xseq, rep(0, length(xseq))))
      } else {
          if(is.character(glue))
            this.cc <- rintcal::glue.ccurves(prebomb=cc, postbomb=glue, cc.dir, as.F=is.F, as.pMC=is.pMC) else
              if(glue>0) {
                if(glue %in% 1:3)
                  this.cc <- rintcal::glue.ccurves(prebomb=1, postbomb=glue, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                    if(glue %in% 4:5)
                      this.cc <- rintcal::glue.ccurves(prebomb=3, postbomb=glue, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                        stop("please provide an integer for glue between 0 and 5")
              } else {
                  this.cc <- rintcal::ccurve(cc, postbomb=postbomb, cc.dir=cc.dir,
                    resample=cc.resample, as.F=is.F, as.pMC=is.pMC)

            # check if any dates lie at the younger edge of this.cc
            young <- FALSE
            if(is.F || is.pMC) {
              if(max(y+(3*er)) > max(this.cc[,2]))
                young <- TRUE
            } else
                if(min(y-(3*er)) < min(this.cc[,2]))
                  young <- TRUE

              if(young) {
                if(postbomb == FALSE) {
                  if(bombalert)
                    stop("please provide a postbomb curve") else {
                      # no postbomb has been defined, but we need to add one
                      if(cc==1) # then assume we need NH1
                        this.cc <- rintcal::glue.ccurves(prebomb=cc, postbomb=1, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                          if(cc==3) # then assume we need SH1-2
                            this.cc <- rintcal::glue.ccurves(prebomb=cc, postbomb=4, cc.dir, as.F=is.F, as.pMC=is.pMC) else
                              stop("please provide a postbomb curve, e.g. postbomb=1")
                    }
                } else
                  this.cc <- rintcal::glue.ccurves(prebomb=cc, postbomb=postbomb, cc.dir, as.F=is.F, as.pMC=is.pMC)
              }
            }
           return(this.cc)
         }
    }
}


