## ----eval=FALSE---------------------------------------------------------------
# install.packages('rice')

## ----eval=FALSE---------------------------------------------------------------
# update.packages()

## -----------------------------------------------------------------------------
library(rice)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
howmuchC14(0)
howmuchC14(55e3)
x <- seq(0, 55e3, length=100) # a sequence of ages
y <- sapply(x, function(x) howmuchC14(x, talk=FALSE))
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
abline(h=BCADtoC14(1694), lty=3)
abline(v=C14toBCAD(129), lty=3)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 2020, BCAD=TRUE, cc2='nh1')

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 2020, BCAD=TRUE, cc2='nh1', add.yaxis=TRUE)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(50000, 35000, timescale="D")

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(cc1="stuiver_suess_1966", cc2="pearson_stuiver_1986")

## -----------------------------------------------------------------------------
calBPtoC14(10.5)
BCADtoC14(1940:1950)

## -----------------------------------------------------------------------------
BCADtocalBP(2025)
BCADtocalBP(-1, zero=TRUE)
BCADtocalBP(-1, zero=FALSE)

## -----------------------------------------------------------------------------
Delta14CtoF14C(152, t=4000)
F14CtoDelta14C(0.71, t=4000)

## -----------------------------------------------------------------------------
C14toDelta14C(152, t=4000)
Delta14CtoC14(592, t=4000)

## ----fig.width=5, fig.asp=1---------------------------------------------------
draw.CF(50000, 3000) 

## ----fig.width=8, fig.asp=.6--------------------------------------------------
fromto(100, "calBP")

## ----fig.width=4, fig.asp=.8--------------------------------------------------
x <- seq(0, 55e3, length=1e3)
cc <- calBPtoC14(x)
Dcc <- C14toDelta14C(cc[,1], cc[,2], x)

par(mar=c(4,3,1,3), bty="l")
plot(x/1e3, Dcc[,1]+Dcc[,2], type="l", xlab="kcal BP", ylab="")
mtext(expression(paste(Delta, ""^{14}, "C")), 2, 1.7)
lines(x/1e3, Dcc[,1]-Dcc[,2])

## ----fig.width=6, fig.asp=.8--------------------------------------------------
draw.ccurve(cc2="IntCal20", timescale2="D", add.yaxis=TRUE)

## -----------------------------------------------------------------------------
data(shroud)
shroud

## -----------------------------------------------------------------------------
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
y <- c(3820, 4430+90)
er <- c(40, sqrt(40^2 + 25^2)) 
cc <- c(1,2)
overlap(y, er, cc=cc)

