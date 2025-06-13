## ----eval=FALSE---------------------------------------------------------------
# install.packages('rice')

## ----eval=FALSE---------------------------------------------------------------
# update.packages()

## -----------------------------------------------------------------------------
library(rice)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
howmanyC14(0)
howmanyC14(55e3)
x <- seq(0, 55e3, length=100) # a sequence of ages
y <- sapply(x, function(x) howmanyC14(x, talk=FALSE))
plot(x, y, type="l", xlab="time (cal BP)", ylab="C-14 remaining")

## -----------------------------------------------------------------------------
adjust.fractionation(9000, -17) # the sample's d13C was measured at -17
adjust.background(9000, 50, 48000, 200) # a well-oiled machine, with backgrounds of F=0.00254

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve()
grid()
abline(0, 1, lty=2)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1000, 2020, BCAD=TRUE, cc2='marine20', add.yaxis=TRUE)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 1950, BCAD=TRUE)

## -----------------------------------------------------------------------------
BCADtoC14(1694) # gives 129 +- 9 14C BP 
C14toBCAD(129) # the curve crosses 129 14C BP 5 times

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 1950, BCAD=TRUE)
abline(h=BCADtoC14(1694)[1], lty=3) 
abline(v=C14toBCAD(129), lty=3)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 2020, BCAD=TRUE, cc2='nh1')

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 2020, BCAD=TRUE, cc2='nh1', add.yaxis=TRUE)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(50000, 35000, realm="D")

## -----------------------------------------------------------------------------
calBPtoC14(10.5)
BCADtoC14(1940:1950)

## -----------------------------------------------------------------------------
BCADtocalBP(2025)
BCADtocalBP(-1, zero=TRUE)
BCADtocalBP(-1, zero=FALSE)

## -----------------------------------------------------------------------------
D14CtoF14C(152, t=4000)
F14CtoD14C(0.71, t=4000)

## -----------------------------------------------------------------------------
C14toD14C(152, t=4000)
D14CtoC14(592, t=4000)

## ----fig.width=8, fig.asp=.6--------------------------------------------------
fromto(100, "calBP")

## ----fig.width=4, fig.asp=.8--------------------------------------------------
x <- seq(0, 55e3, length=1e3)
cc <- calBPtoC14(x)
Dcc <- C14toD14C(cc[,1], cc[,2], x)

par(mar=c(4,3,1,3), bty="l")
plot(x/1e3, Dcc[,1]+Dcc[,2], type="l", xlab="kcal BP", ylab="")
mtext(expression(paste(Delta, ""^{14}, "C")), 2, 1.7)
lines(x/1e3, Dcc[,1]-Dcc[,2])

## ----fig.width=6, fig.asp=.8--------------------------------------------------
draw.ccurve(cc2="IntCal20", realm2="D", add.yaxis=TRUE)

## -----------------------------------------------------------------------------
data(shroud)
shroud
pool(shroud$y,shroud$er) 
Zu <- grep("ETH", shroud$ID) # Zurich lab only
pool(shroud$y[Zu],shroud$er[Zu])

## -----------------------------------------------------------------------------
weighted_means(shroud$y[Zu],shroud$er[Zu])

## ----fig.width=5, fig.asp=.8--------------------------------------------------
as.one(shroud$y,shroud$er)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
as.bin(shroud$y,shroud$er, 100, move.by=25)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
spread(shroud$y,shroud$er)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
span(700, 20, 750, 20)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
y <- c(3820, 4590-90) 
er <- c(40, sqrt(40^2 + 25^2)) 
cc <- c(1,2)
overlap(y, er, cc=cc)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
contaminate(5000, 20, 10, 0, 1)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
contaminate(66e6, 1e6, 0.5, 0.1, 1)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
clean(9000, 100, 10)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
muck(591, 30, BCADtoC14(40)[,1], 0, 1)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
perFaith <- BCADtoC14(40)
repairF <- BCADtoF14C(1400)
muck(591, 30, perFaith[,1], perFaith[,2], repairF[,1], repairF[,2])

## ----fig.width=6, fig.asp=.8--------------------------------------------------
real.14C <- seq(0, 50e3, length=200)
contam <- seq(0, 10, length=101) # 0 to 10% contamination
contam.col <- rainbow(length(contam))
plot(0, type="n", xlim=c(0, 55e3), xlab="real 14C age", ylim=range(real.14C), ylab="observed 14C age")
for (i in 1:length(contam)) {
  observed <- contaminate(real.14C, 0, contam[i], 0, 1, talk=FALSE, MC=FALSE)
  lines(real.14C, observed[,1], col = contam.col[i])
}

## ----fig.width=6, fig.asp=.8--------------------------------------------------
draw.contamination()

## -----------------------------------------------------------------------------
Cs <- c(.02, .05, .03, .04) # carbon contents of each fraction
wghts <- c(5, 4, 2, .5) # weights for all fractions, e.g., in mg
ages <- c(130, 130, 130, NA) # ages of all fractions. The unmeasured one is NA
errors <- c(10, 12, 10, NA) # errors, unmeasured is NA
fractions(150, 20, Cs, wghts, ages, errors) # assuming a bulk age of 150 +- 20 C14 BP 

## ----fig.width=4, fig.asp=.8--------------------------------------------------
calib.130 <- caldist(130, 10, BCAD=TRUE)
plot(calib.130, type="l")

## -----------------------------------------------------------------------------
l.calib(145, 130, 10)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
dice <- r.calib(100, 130, 10)
plot(density(dice))
rug(dice)

## -----------------------------------------------------------------------------
hpd(calib.130)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
calib.2450 <- caldist(2450, 20)
plot(calib.2450, type="l")
points.2450 <- point.estimates(calib.2450)
points.2450
abline(v=points.2450, col=1:4, lty=2)

## ----fig.width=5, fig.asp=1---------------------------------------------------
calibrate(2450, 40)

## ----fig.width=5, fig.asp=1---------------------------------------------------
mycurve <- smooth.ccurve(smooth=50)
calibrate(2450, 40, thiscurve=mycurve)

## ----fig.width=5, fig.asp=1---------------------------------------------------
try(calibrate(130,30))
calibrate(130, 30, bombalert=FALSE)

## -----------------------------------------------------------------------------
younger(150, 130, 10)
older(150, 130, 10)

## -----------------------------------------------------------------------------
p.range(300, 150, 130, 10)

## ----fig.width=5, fig.asp=1---------------------------------------------------
push.normal(2450,40, 400,20)

