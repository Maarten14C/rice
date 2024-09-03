# rice 0.1.2
* added an option bombalert to the calibrate function. If set to false, it can plot ages close to 0 C14 BP without warnings.
* added the data from the marine database (calib.org/marine), as data 'shells'
* new functions find.shells and map.shells to plot shells data in maps based on their coordinates
* new function shells.mean to plot deltaRs of selected shells, and calculate a weighted mean deltaR
* new function weighted_means to calculate weighted means and errors for multiple radiocarbon dates (or delta R values)

# rice 0.1.1
* added citation information
* added a function 'older'
* added a vignette

# rice 0.1.0
* This is the first release of the rice package. It now separates the calibration functions from its parent data package rintcal, which in the future will contain the IntCal and other calibration curve data only.
