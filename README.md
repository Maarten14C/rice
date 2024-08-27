# rice

## Radiocarbon Calibration Equations

Radiocarbon is widely used for dating a range of archaeological and geographical objects. This package provides a number of functions to calibrate radiocarbon dates, including ones to transfer values between different realms (F14, D14C, pMC, C14 ages), to simulate the impacts of contamination, and to plot one or more calibrated dates.

Please check out the vignettes folder for a tutorial. In short, install the package from github, load it and run the two main functions:

```{r, eval=FALSE}
require(devtools)
install_github("Maarten14C/rice")
library(rice)
calibrate(2450,20)
contaminate(55000, 100, .05, 1)
draw.contamination()
```
