## ----eval=FALSE---------------------------------------------------------------
#  install.packages('rice')

## ----eval=FALSE---------------------------------------------------------------
#  update.packages()

## -----------------------------------------------------------------------------
library(rice)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve()

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1000, 2020, BCAD=TRUE, cc2='marine20', add.yaxis=TRUE)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
draw.ccurve(1600, 1950, BCAD=TRUE)

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
BCADtocalBP(2024)
BCADtocalBP(-1, zero=TRUE)
BCADtocalBP(-1, zero=FALSE)

## -----------------------------------------------------------------------------
D14CtoF14C(152, t=4000)
F14CtoD14C(0.71, t=4000)

## -----------------------------------------------------------------------------
C14toD14C(0.71, t=4000)
D14CtoC14(152, t=4000)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
x <- seq(0, 55e3, length=1e3)
cc <- calBPtoC14(x)
Dcc <- calBPtoD14C(x)

par(mar=c(4,3,1,3), bty="l")
plot(x/1e3, Dcc[,1]+Dcc[,2], type="l", xlab="kcal BP", ylab="")
mtext(expression(paste(Delta, ""^{14}, "C")), 2, 1.7)
lines(x/1e3, Dcc[,1]-Dcc[,2])

par(new=TRUE)
plot(x/1e3, (cc[,1]-cc[,2])/1e3, type="l", xaxt="n", yaxt="n", col=4, xlab="", ylab="")
lines(x/1e3, (cc[,1]+cc[,2])/1e3, col=4)
mtext(expression(paste(""^{14}, "C kBP")), 4, 2, col=4)
axis(4, col=4, col.axis=4, col.ticks=4)

## -----------------------------------------------------------------------------
data(shroud)
shroud
pool(shroud$y,shroud$er) 
Zu <- grep("ETH", shroud$ID) # Zurich lab only
pool(shroud$y[Zu],shroud$er[Zu])

## ----fig.width=5, fig.asp=.8--------------------------------------------------
as.one(shroud$y,shroud$er)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
as.bin(shroud$y,shroud$er, 50, 10)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
spread(shroud$y,shroud$er)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
contaminate(5000, 20, 1, 1)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
contaminate(66e6, 1e6, 0.5)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
clean(9000, 100, percentage=1)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
muck(591, BCADtoC14(40)[1], 1)

## ----fig.width=5, fig.asp=.8--------------------------------------------------
muck(591, BCADtoC14(40)[1], BCADtoF14C(1400)[1])

## ----fig.width=6, fig.asp=.8--------------------------------------------------
real.14C <- seq(0, 50e3, length=200)
contam <- seq(0, 10, length=101) # 0 to 10% contamination
contam.col <- rainbow(length(contam))
plot(0, type="n", xlim=c(0, 55e3), xlab="real 14C age", ylim=range(real.14C), ylab="observed 14C age")
for (i in 1:length(contam)) {
  observed <- contaminate(real.14C, 0, contam[i], 1, decimals=5, talk=FALSE)
  lines(real.14C, observed[,1], col = contam.col[i])
}
contam.legend <- seq(0, 10, length=6)
contam.col <- rainbow(length(contam.legend)-1)
text(50e3, contaminate(50e3, 0, contam.legend, 1, visualise=FALSE, talk=FALSE)[,1],
  labels=contam.legend, cex=.7, offset=0, adj=c(0,.8))

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

## -----------------------------------------------------------------------------
hpd(calib.130)

## ----fig.width=4, fig.asp=.8--------------------------------------------------
calib.2450 <- caldist(2450, 20)
plot(calib.2450, type="l")
points.2450 <- point.estimates(calib.2450)
points.2450
abline(v=points.2450, col=1:4, lty=2)

## ----fig.width=5, fig.asp=1---------------------------------------------------
calibrate(2450, 20)

## ----fig.width=5, fig.asp=1---------------------------------------------------
mycurve <- smooth.ccurve(smooth=50)
calibrate(2450, 20, thiscurve=mycurve)

## ----fig.width=5, fig.asp=1---------------------------------------------------
try(calibrate(130,30))
calibrate(130, 30, bombalert=FALSE)

## -----------------------------------------------------------------------------
younger(150, 130, 10)
older(150, 130, 10)

## ----fig.width=5, fig.asp=1---------------------------------------------------
myshells <- map.shells(S=54, W=-8, N=61, E=0) # the northern part of the UK

## ----fig.width=5, fig.asp=1---------------------------------------------------
head(myshells)
shells.mean(myshells)

## ----fig.width=5, fig.asp=1---------------------------------------------------
myshells <- find.shells(120, 10, 20)
shells.mean(myshells, distance=TRUE)

## ----fig.width=4, fig.asp=1---------------------------------------------------
set.seed(123)
dates <- sort(sample(500:2500,5))
errors <- .05*dates
depths <- 1:length(dates)
my.labels <- c("my", "very", "own", "simulated", "dates")
draw.dates(dates, errors, depths, BCAD=TRUE, labels=my.labels, age.lim=c(0, 1800))

## ----fig.width=4, fig.asp=1---------------------------------------------------
plot(300*1:5, 5:1, xlim=c(0, 1800), ylim=c(5,0), xlab="AD", ylab="dates")
draw.dates(dates, errors, depths, BCAD=TRUE, add=TRUE, labels=my.labels, mirror=FALSE)

## ----fig.width=4, fig.asp=1---------------------------------------------------
par(bg="black", mar=rep(1, 4))
n <- 50; set.seed(1)
draw.dates(rnorm(n, 2450, 30), rep(25, n), n:1,
  mirror=FALSE, draw.base=FALSE, draw.hpd=FALSE, col="white",
  threshold=1e-28, age.lim=c(2250, 2800), ex=.8)

