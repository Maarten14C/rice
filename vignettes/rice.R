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

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve()

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

## ----fig.width=4, fig.asp=.8--------------------------------------------------
fromto(100, "calBP")

## ----fig.width=4, fig.asp=.8--------------------------------------------------
x <- seq(0, 55e3, length=1e3)
cc <- calBPtoC14(x)
Dcc <- C14toD14C(cc[,1], cc[,2], x)

par(mar=c(4,3,1,3), bty="l")
plot(x/1e3, Dcc[,1]+Dcc[,2], type="l", xlab="kcal BP", ylab="")
mtext(expression(paste(Delta, ""^{14}, "C")), 2, 1.7)
lines(x/1e3, Dcc[,1]-Dcc[,2])

## ----fig.width=4, fig.asp=.8--------------------------------------------------
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

