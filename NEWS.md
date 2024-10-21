# rice 0.3.0
* hpd ranges are now calculated at a specified precision (defaults to yearly)
* new function `BCADto14C` to calculate 14C age belonging to a BC/AD age (this calls the function `calBPto14C`)
* new option `as.F` to calibrate in the F14C realm (the default remains to calibrate in the C14 realm)
* warnings are printed on calibrate() plots if dates are truncated and edge=TRUE (with the default edge=FALSE, dates that are truncated are not calibrated). The printed warning can be removed by setting `print.truncate.warning=FALSE`
* renaming `age.F14C`, `F14C.age`, `age.pMC` and `pMC.age` to, respectively, `C14toF14C`, `F14CtoC14`, `C14topMC` and `pMCtoC14`. This because `age` is an ambivalent term in this context
* new functions `BCADtocalBP` and `calBPtoBCAD` to transfer cal BP into BC/AD ages and vice versa. Can deal with (e.g. Gregorian/Julian) calendars which do not include 0
* new functions to translate between any of the realms `calBP`, `BCAD`, `C14`, `F14C`, `pMC` and `D14C`.
* new function `smooth.ccurve` to smooth a calibration curve using a moving window of a specified width. This can be useful to calibrate material that is known to have accumulated, say, over two decades
* new function `pool` which calculates the chi2 and accompanying p-value for a set of multiple measurements on the same sample. If the scatter between the values is low enough for the p-value to be below a threshold, then the pooled mean and uncertainty are returned
* the function draw.dates now has an option `oncurve` to draw the dates onto the calibration curve.
* added dataset `shroud`, which contains replicate radiocarbon measurements on the Shroud of Turin, from three labs.
* new function `decontaminate` to estimate the percentage of contamination needed to explain the difference between a 'real' and an 'observed' radiocarbon age.

# rice 0.2.0
* added an option `bombalert` to the calibrate function. If set to false, plots ages close to 0 C14 BP without warnings.
* added the data from the marine database (calib.org/marine), as data `shells`
* new functions `find.shells` and `map.shells` to plot shells data in maps based on their coordinates
* new function `shells.mean` to plot deltaRs of selected shells, and calculate a weighted mean deltaR
* new function `weighted_means` to calculate weighted means and errors for multiple radiocarbon dates (or delta R values)
* repaired a bug in `draw.D14C`
* `draw.ccurve` now can plot the C14 in the 'realms' of C14 BP, F14C, pMC and D14C using the 'realm' option.

# rice 0.1.1
* added citation information
* added a function `older`
* added a vignette

# rice 0.1.0
* The first release of the rice package. It separates the calibration functions from its parent data package rintcal, which in the future will contain the IntCal and other calibration curve data only.
