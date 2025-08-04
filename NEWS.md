# rice 1.4.0
* renamed 'realm' to the more accurate term 'timescale'.
* `BCADtoD14C` and `b2ktoD14C` now calculate D14C values correctly.

# rice 1.3.0
* `caldist` now deals better with uncalibrated distributions (`cc=0`).
* Where relevant, timescale-related functions such as `C14toF14C` now have more flexibility regarding rounding of values. The default is no rounding. 
* new function `draw.CF` to visualise the relationship between 14C ages and their F values, and how ages with large errors become skewed/asymmetric.
* the functions which use random samples (e.g., `contaminate`) now have the option for a seed to be set.
* the function `muck` now not only returns the percentage of contamination required to to from an observed age to a target one (and an assumed F14C of the contamination). As an alternative, if the percentage of contamination is known, then the required F14C of the contamination can be returned.
* The functions `F14CtoC14` and `C14toF14C` can now also be called with shorter names - `FtoC` and `CtoF`, respectively.
* `draw.ccurve` now plots the title of the second axis (when `add.yaxis=TRUE`).
* new function `coverage` which calculates to what degree a distribution A is covered by (i.e., falls within) another distribution B. This could be used to check for example how well an age estimate from an age-model fits with a calibrated date. If both overlap, even if the date has a much wider distribution than the model estimate, then the coverage is high (even if the total overlap between the narrow and wider distributions is much less than 100\%). 
* new function `hpd.overlap` which checks if any of the highest posterior density (hpd) intervals of two distributions overlap. Returns TRUE if any of the hpd ranges overlap, FALSE if not.
* the function `fractions` now calculates a combined age if the weights, percentage carbon and ages of all fractions are provided. 

# rice 1.2.0
* new timescale functions to convert from/to b2k (years before AD 2000), an age scale popular in the ice core community.
* new optional browsable plots for `find.shells` and `map.shells` functions. Use `browse=TRUE`. Requires Internet connection. 
* The `find.shells` and `map.shells` functions now also offer to plot a browsable maps of (current) ocean currents. This to help interpret which shells would be most representative of different ocean water masses. Requires Internet connection. 
* titles of contaminate plots are now plotted inside of the device range.
* new functions `adjust.fractionation` and `adjust.background` to correct for fractionation and background values (still experimental).
* the `calibrate` function now draws correctly when is.pMC=TRUE.
* the `overlapping` function can now calculate the overlap between either radiocarbon ages or distributions (the latter should be provided as lists).
* the function `D14CtoC14` had a bug that misinterpreted columns - should work better now.
* added testthat functions to facilitate debugging.

# rice 1.1.1
* added a secondary y axis to the `contaminate`, `clean` and `muck` plots.
* additional R packages for plotting maps are now only installed and loaded when required.
* added the latest entries from the calib.org/marine database (2 April 2025). 
* draw.dates now doesn't throw an error when ka=TRUE, and also works better when rotate.axes=TRUE.
* hpd gains an option `bins` to provide the minimum number of bins in a distribution before hpds are calculated. Any distributions with fewer bins get recalculated using a narrower binsize (equating to 100 equally-spaced bins).
* added a function `span` to calculate (calibrated) time-spans between two radiocarbon dates.
* repaired a bug in C14tocalBP which caused NAs for some cal BP values.
* a function `overlapping` calculates to what degree two calibrated dates are overlapping.
* within the timescale-related functions, changed the rounding procedure from `signif` to `round`. 
* added a function `fromto` that translates values into different domains, and plots them.
* calibrated distributions are now drawn in more consistent ways (e.g., direction up/down, height).

# rice 1.0.0
* corrected a bug in `C14topMC` where errors were not calculated correctly.
* caldist now always glues a postbomb curve to a prebomb if postbomb is not FALSE.
* new functions `push.normal` and `push.gamma` to push a date to younger or older ages by adding/subtracting a normal resp. gamma distribution.
* reservoir effect (deltaR, deltaSTD) options have been added to functions where this is relevant.
* added more warnings to `calibrate` for very old dates.
* calibrated distributions are now plotted more consistently between functions.
* added uncertainty estimates to the `clean`, `contaminate` and `muck` functions, using Monte Carlo-based sampling.
* Some of the timescale functions translating to D14C space didn't handle multiple entries well. This should work correctly now.
* `find.shells` and `map.shells` now deal better with missing mapping-related packages.

# rice 0.4.0
* the `contaminate` function now also produces a plot, and more details of the calculations.
* the `decontaminate` function has been renamed to `clean` and has been updated with clearer messages; it now removes contamination to calculate the true/target age. A plot is made, and calculation details are provided.
* new function `muck` to calculate how much contamination has to be inferred to go from an observed age to a true/target age. A plot is made, and calculation details are provided.
* new functions `C14tocalBP` and `C14toBCAD`. Since these rely on the outdated 'intercept calibration' method, these functions are provided for illustrative purposes only.
* draw.dates can now also plot the dates on the calibration curve, using the `oncurve` option. If so, then the curve and dates can also be plotted in the `F14C` or `pMC` timescales.
* New function `r.calib` that samples random calendar ages from a calibrated distribution.
* new function `p.range` to calculate a calibrated age's probability of lying between a range of BC/AD or cal BP ages.
* new function `as.one` to calculate the product of multiple calibrated ages, assuming that they all stem from (exactly) the same calendar age. Not that this is dangerous, and care should be taken to make sure that the assumptions are met.
* new function `as.bin` to calculate how many of a set of calibrated radiocarbon dates fall into bins of a specified width. The bin moves along the range of calibrated ages, to visualise how many dates fit bins over time. This would be safer than using the function `as.one`. 
* new function `spread` shows the spread (in calendar years) of a set of dates. Accompanies the functions `pool` and `as.one`.

# rice 0.3.0
* hpd ranges are now calculated at a specified precision (defaults to yearly)
* new function `BCADto14C` to calculate the 14C age belonging to a BC/AD age (this calls the function `calBPto14C`)
* new option `as.F` to calibrate in the F14C timescale (the default remains to calibrate in the C14 timescale)
* warnings are printed on calibrate() plots if dates are truncated and edge=TRUE (with the default edge=FALSE, dates that are truncated are not calibrated). The printed warning can be removed by setting `print.truncate.warning=FALSE`
* renaming `age.F14C`, `F14C.age`, `age.pMC` and `pMC.age` to, respectively, `C14toF14C`, `F14CtoC14`, `C14topMC` and `pMCtoC14`. This because `age` is an ambivalent term in this context
* new functions `BCADtocalBP` and `calBPtoBCAD` to transfer cal BP into BC/AD ages and vice versa. Can deal with (e.g. Gregorian/Julian) calendars which do not include 0
* new functions to translate between any of the timescales `calBP`, `BCAD`, `C14`, `F14C`, `pMC` and `D14C`.
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
* `draw.ccurve` now can plot the C14 in the timescales of C14 BP, F14C, pMC and D14C using the `timescale` option.

# rice 0.1.1
* added citation information
* added a function `older`
* added a vignette

# rice 0.1.0
* The first release of the rice package. It separates the calibration functions from its parent data package rintcal, which in the future will contain the IntCal and other calibration curve data only.
