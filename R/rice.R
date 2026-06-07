
# can intcal.data be made to show the SH data? Currently shows the NH data

# sample weight functions (per Philippa Ascough's suggestion). Given a %C (perhaps provide estimates for sample types such as peat, bone, ...), a loss during pretreatment, and a required graphite weight, what sample weight will be required?)

# prepare a function to redo deltaR calcs when new Marine curves come out. Using BCADtocalBP(shells$collected), calBPto14C(cc=2) and shells$C14, shells$er. Unclear how the dR errors are obtained.

# fruits-type model that mixes atmospheric and marine calibration curves. Freshwater effects can cause C14 shifts of up to 1k.

# error multipliers, rounding. Could add procedures for different labs, e.g. QUB_bg, etc. This would be useful for reasons of transparency and community standards. Add data from historical UBA standards/backgrounds?



#' @name howmuchC14
#' @title Amount of C14 particles in a sample
#' @description Calculate the expected amount of remaining C14 atoms in a sample, given its weight and age.
#' @details The number of carbon atoms in the sample is estimated. Given the known C14/C ratio at F=1, and given the sample's age, we can estimate the number of remaining C14 atoms. Given a 12C current at the detector end of an AMS, we can then also calculate how many 14C ions would be counted per second and minute. 
#' Note that backgrounds are not modelled (but could be investigated by e.g. typing \code{howmuchC14(45e3)} which gives as c. 1 background count per second). The calculated C14 count rate assumes no fractionation.
#' @return The estimated number of C14 atoms.
#' @param age The age of the sample (in cal BP per default, or in C14 BP if use.cc=FALSE).
#' @param wght The weight of the sample (in mg). Defaults to 1 mg.
#' @param use.cc Whether or not to use the calibration curve. If set to \code{use.cc=FALSE}, then we assume that the age is the radiocarbon age (this enables ages beyond the reach of the calibration curves to be used).
#' @param Av Avogadro's number, used to calculate the number of carbon atoms in the sample.
#' @param C14.1950 The standard 14C/C ratio back in AD 1950 (1.176e-12, so around 1 in 1 trillion carbon atoms was a 14C atom at that moment in time.
#' @param current The current of 12C+ ions arriving at the Faraday counter. Defaults to \code{current=25e-6}, 25 micro-Ampere.
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
howmuchC14 <- function(age, wght=1, use.cc=TRUE, Av=6.02214076e23, C14.1950=1.176e-12, current=25e-6, format="g", cc=1, postbomb=FALSE, glue=0, cc.dir=NULL, thiscurve=NULL, talk=TRUE, as.AMS=TRUE, decimals=3) {
  if(length(age) > 1)
    stop("can only handle single value for age")
  
  if(use.cc) {
    F <- calBPtoF14C(age, cc=cc, postbomb=postbomb, glue=glue, cc.dir=cc.dir, thiscurve=thiscurve)[,1]
    if(is.na(F)) {
      message("Cannot use calibration curve for this age, assuming C14 age")
      F <- C14toF14C(age)
    }
  } else
      F <- C14toF14C(age) # then 'age' is on the C14 scale
  F <- as.numeric(F)
  
  atoms <- (wght/1e3)*Av/12 # number of C atoms in a mg
  C14_atoms <- F * C14.1950 * atoms # number of C14 atoms in 1950 in the same amount of sample
  C14 <- as.numeric(ceiling(F * C14.1950 * atoms)) # C14 atoms roundest up
  
  # from the current, calculate amount of ions ending up at the C12 Faraday cup:
  i12 <- current/1.602e-19 # e, electric charge for a single proton, in coulomb
  # from this, we can model the numbers of C14 atoms arriving at the C14 detector:
  i14 <- C14.1950 * i12 * F # charge in A = protons arriving per second
  persecond <- round(i14, decimals)
  i12 <- formatC(i12, format=format, digits=decimals)
  
  atoms <- formatC(atoms, format=format, digits=decimals)
  C14.talk <- formatC(C14, format=format, digits=decimals)
  decays <- formatC(C14 * log(2) / (5730 * 365.25), format="fg", digits=decimals)

  if(talk) {
    message(wght, " mg carbon contains ", atoms, " C atoms")
    message("14C atoms remaining at ", age, " cal BP (F=", round(F, decimals), "): ", C14.talk)
    message(decays, " 14C atoms in the sample will decay each day")
	if(as.AMS)
      message(paste0("For a 12C current of ", current*1e6, " micro-ampere (", i12,
        " 12C/second) at the AMS detector,\n  ",
        persecond, " 14C particles would be counted per second (",
        60*persecond, " per minute)" ))
  }

  invisible(list(C14=C14, C14persecond=i14, decaysperday=as.numeric(decays), i12=current/1.602e-19, C14.1950=C14.1950, totC=as.numeric(atoms)))
}



#' @name radio
#' @title Listen to radiocarbon being measured
#' @description A sonification of the amount of C14 being detected in an AMS. 
#' @details The sonification can be heard as 'Geiger counter-like' clicks (individual C14 atoms ending up in the detector) and/or as a possibly slightly meandering sine wave (e.g., for a sample of age 5000 14C BP, we'd count c. 105.6 atoms per second in the AMS, so this would become a tone of frequency 105.6 Hz). The calculations are based on the function \code{howmuchC14}.
#'
#' Instead of counting the particles themselves (as done in an AMS), we can also count the decay events as used in radiometric counting. In this case, we need much larger samples (wght) and longer counting times (duration).
#'
#' Note that the calculations in this function do not include error multipliers or corrections such as for fractionation and backgrounds. 
#' @return A sound (tone and or clicks) is played and returned invisibly.
#' @param age The age of the sample (in cal BP per default, or in C14 BP if use.cc=FALSE).
#' @param duration How long the sample will count (and sound) for in seconds. Defaults to 10 seconds.
#' @param as.decays Work with the C14 decays. Defaults to FALSE, which works with the number of C14 atoms instead. If set to true, you'll need larger sample sizes (wght) and longer counting times (duration) to get decent counts.
#' @param wght The weight of the sample (in mg). Defaults to 1 mg.
#' @param as.clicks Make the C-14 counts sound as clicks, based on random Poisson sampling.
#' @param as.tone Make the C-14 frequency (counts per second) sound as a wave (e.g., 105 cps becomes a 105 Hz sine wave). This can either be a constant wave or be meandering (see `noise` and `meander`).
#' @param click_length Length of the clicks. Defaults to 80, but can be altered to make the clicks sound different
#' @param tone.volume Volume of the tone/wave relative to that of the clicks.
#' @param noise Deviation from the mean of the tone. Defaults to 0.02.
#' @param meander Drift of the tone along the mean. Defaults to 0.005. 
#' @param play Whether or not to play the sound
#' @param sr Sampling rate. This audio quality defaults to 44100 (per second), which is based on the CD standard.
#' @param visualise the counts and calculation of the C14 age. Defaults to TRUE.
#' @param cex Size of the font. Defaults to \code{cex=0.6}.
#' @param return.sound Return the sound as an invisible object. If set to FALSE (the default), instead returns the counts and the calculated C14 age.
#' @param ... Optional constants to be entered into the function `howmuchC14`
#' @author Maarten Blaauw
#' @examples
#'   radio(0)
#'   radio(45000)
#'   # decay events over 1 minute in 1 gram of carbon of age 500 14C BP:
#'   radio(500, wght=1000, as.decays=TRUE, duration=60) 
#' @export
radio <- function(age, duration=10, as.decays=FALSE, wght=1, as.clicks=TRUE, as.tone=TRUE, click_length=80, tone.volume=0.5, noise=0.02, meander=0.005, play=interactive(), sr=44100, visualise=TRUE, cex=.6, return.sound=FALSE, ...) {

  hasaudio <- requireNamespace("audio", quietly=TRUE)
  hastuneR <- requireNamespace("tuneR", quietly=TRUE)
  if(!hasaudio)
    stop("Please install the audio package:\ninstall.packages(\"audio\")")

  if(as.decays)
	as.tone <- FALSE # no music

  counts <- howmuchC14(age, wght=wght, as.AMS=!as.decays, ...)
  if(as.decays) # rate = per second
    rate <- counts$decaysperday / (60*60*24) else # C14 decays
      rate <- counts$C14persecond # C14 atoms
	  
  n <- sr * duration # number of steps
  p <- rate / sr # probability per sample
  events <- runif(n) < p # counts/decays, simulated by random Bernoulli trials

  if(visualise) {
	plot(0, type="n", xlim=c(0, duration), ylim=c(0, 1), xlab="duration (s)", ylab="", yaxt="n", yaxs="i", bty="n")
	alpha <- 1/(1+((rate*duration)/200)^0.95) # at higher rates, lines become more transparent
	abline(v=which(events)/sr, col=rgb(0,0,0, alpha))
	rect(-10, .5, duration+10, 1.1, col="white", border=NA)
    # C14_count <- rate*duration
	C14_count <- sum(events)
	
    if(as.decays) { # radiometric: counting decays, not atoms
      C12_count <- counts$totC
      AN <- counts$C14 / counts$totC
	  asF <- AN / counts$C14.1950
	  if(C14_count > 0) {
        AN.sd <- AN / sqrt(C14_count)
        asF.sd <- asF / sqrt(C14_count)
        as.C14 <- F14CtoC14(asF, asF.sd, roundby=0)
      } else {
          AN.sd <- NA
          asF.sd <- NA
          as.C14 <- F14CtoC14(asF, roundby=0)
        }
    } else { # AMS: counting atoms
        C12_count <- counts$i12 * duration
        AN <- C14_count / C12_count
        asF <- AN / counts$C14.1950
        if(C14_count > 0) {
          AN.sd <- AN * sqrt(1/C14_count + 1/C12_count)
          asF.sd <- AN.sd / counts$C14.1950
		  as.C14 <- F14CtoC14(asF, asF.sd, roundby=0)
        } else {
            AN.sd <- NA
            asF.sd <- NA
			as.C14 <- F14CtoC14(asF, roundby=0)
          }
      }
	
    if(as.decays)
      labels <- c(expression(C14~decays~(A %+-% sqrt(A))),
        expression(C12~amount~(N %+-% sqrt(N))), 
        "C14/C12 ratio (AN):", "relative to AD 1950 (F=AN/A0):", 
        "as C14 BP (-8033*log(F)):") else
          labels <- c(expression(C14~counts~(N %+-% sqrt(N))),
            expression(C12~counts~(N %+-% sqrt(N))), 
              "C14/C12 ratio (AN):", "relative to AD 1950 (F=AN/A0):", 
              "as C14 BP (-8033*log(F)):")
	values <- c(
      paste(format(round(C14_count), format="fg"), "\u00B1", format(round(sqrt(C14_count)), format="fg")), 
      paste(format(round(C12_count), format="fg"), "\u00B1", format(round(sqrt(C12_count)), format="fg")),
      paste(format(AN, scientific=TRUE), "\u00B1", format(AN.sd, scientific=TRUE)),
      paste(format(asF, format="fg"), "\u00B1", format(asF.sd, format="fg")),
      paste(as.C14[1], "\u00B1", as.C14[2]))
    legend(0, 1, legend=c(labels, values), cex=cex, ncol=2, 
	  text.col=rep(c(1,rep(grey(.5), 3), 1),2), bty="n")
  }
	  
  # make random clicks based on the rate
  click <- rnorm(click_length) * exp(-seq(0, 6, length.out = click_length))
  click <- click / max(abs(click))
  clicks <- numeric(n)
  for(i in which(events)) {
    end <- min(i + click_length - 1, n)
    clicks[i:end] <- clicks[i:end] + click[1:(end - i + 1)]
  }
  if(max(abs(clicks)) > 0)
    clicks <- clicks / max(abs(clicks)) # normalise

  # make a sine wave based on the rate, perhaps wobble it a bit
  tone <- sin(2 * pi * rate * 1:(n) / sr)
  if(noise > 0) {
    noise <- rnorm(n, 0, noise * sqrt(rate))
    noise[is.na(noise)] <- 0
    wobble <- numeric(n)
    wobble[1] <- rate
    for(i in 2:n)
      wobble[i] <- max(0, wobble[i-1] + meander * (rate - wobble[i-1]) + noise[i])
    tone <- sin(cumsum(2 * pi * wobble / sr))
  }

  if(as.clicks && as.tone)
    sound <- rbind(tone.volume*tone, clicks) else
      if(as.clicks)
        sound <- clicks else 
          sound <- tone

  if(play && (as.clicks+as.tone>0))
    if(nrow(audio::audio.drivers()) > 0) # expected drivers are installed
      audio::play(audio::audioSample(sound, sr)) else {
        if(!hastuneR)
          stop("Please install the tuneR package:\ninstall.packages(\"tuneR\")")
        if(is.matrix(sound))
          sound <- sound[1,] + sound[2,]
        sound <- sound / max(abs(sound))
        sound <- as.integer(sound * 32767) # scale to int16 range (full range of 16 bit sound)

        wave <- tuneR::Wave(left=sound, samp.rate=sr, bit=16)
        tmpwav <- tempfile(fileext=".wav")
        tuneR::writeWave(wave, tmpwav)
        system2("ffplay", args=c("-nodisp", "-autoexit", "-loglevel", "quiet", tmpwav))
      }

  if(return.sound)
    invisible(sound) else
      invisible(list(C14_count=C14_count, duration=duration, F14C=asF, C14_age=as.C14))
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


