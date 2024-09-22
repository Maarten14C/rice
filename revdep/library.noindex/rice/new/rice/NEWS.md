# rice 0.1.2
* added an option bombalert to the calibrate function. If set to false, it can plot ages close to 0 C14 BP without warnings.
* added the data from the marine database (calib.org/marine), as data 'shells'
* new functions find.shells and map.shells to plot shells data in maps based on their coordinates
* new function shells.mean to plot deltaRs of selected shells, and calculate a weighted mean deltaR
* new function weighted_means to calculate weighted means and errors for multiple radiocarbon dates (or delta R values)
* repaired a bug in draw.D14C
* draw.ccurve now can plot the C14 in the 'realms' of C14 BP, F14C, pMC and D14C using the 'realm' option.  

# rice 0.1.1
* added citation information
* added a function 'older'
* added a vignette

# rice 0.1.0
* The first release of the rice package. It separates the calibration functions from its parent data package rintcal, which in the future will contain the IntCal and other calibration curve data only.
