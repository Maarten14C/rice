test_that("howmanyC14 returns expected C14 atoms", {
  result <- howmanyC14(55e3, talk=FALSE)
  expect_length(result, 1)
  expect_equal(unname(result), 115446, tolerance = 1e-5)
})

test_that("adjust.fractionation returns expected adjusted C14 age", {
  result <- adjust.fractionation(5000, -17)
  expect_length(result, 1)
  expect_equal(result, 5131.286, tolerance = 1e-5)
})

test_that("adjust.background returns expected background-corrected C14 age", {
  result <- unname(unlist(adjust.background(9000, 50, 45000, 200)))
  expect_length(result, 2)
  expect_equal(result[1], 9061.71822, tolerance = 1e-5)
  expect_equal(result[2], 50.94971, tolerance = 1e-5)
})

### marine

test_that("find.shells returns expected amount of shells", {
  myshells <- find.shells(0, 55, mapsize="small", currents=FALSE)
  expect_length(myshells, 16)
})

test_that("weighted_means returns expected weighted mean of selected shells", {
  myshells <- find.shells(0, 55, mapsize="small", currents=FALSE)
  wmns <- weighted_means(myshells[,5], myshells[,6], talk=FALSE)
  expect_length(wmns, 2)
  expect_equal(wmns, c(-155.2, 66.2), tolerance = 1e-5)
})

test_that("shells.mean returns expected mean of shells", {
  myshells <- find.shells(0, 55, mapsize="small", currents=FALSE)
  wmns <- shells.mean(myshells, talk=FALSE)
  expect_length(wmns, 2)
  expect_equal(wmns, c(-155.2, 66.2), tolerance = 1e-5)
})

### calibrate

test_that("point.estimates returns expected point estimates", {
  result <- unname(point.estimates(caldist(130,20)))

  expect_length(result, 4)
  expect_equal(result, c(128.2, 105.8, 76.0, 139.5), tolerance = 1e-5)
})


test_that("point.estimates returns expected point estimates", {
  result <- hpd(caldist(130,20))
  expect_length(result, 3)
  expect_equal(result[,3], c(68.7, 1.9, 24.5), tolerance = 1e-5)
})


test_that("l.calib returns expected likelihood", {
  result <- l.calib(100, 130, 20)
  expect_length(result, 1)
  expect_equal(result, 0.01721038, tolerance = 1e-5)
})

test_that("r.calib returns expected random points", {
  set.seed(123)
  result <- r.calib(10, 130, 20)
  expect_length(result, 10)
  expect_equal(result[1], 74.31000, tolerance = 1e-5)
})

test_that("older returns expected probability", {
  result <- older(2500, 2450, 20)
  expect_length(result, 1)
  expect_equal(result, 0.5302995, tolerance = 1e-5)
})

test_that("p.range returns expected probability", {
  result <- p.range(2800, 2400, 2450, 20)
  expect_length(result, 1)
  expect_equal(result, 0.9127575, tolerance = 1e-5)
})


test_that("calib.t returns expected ages", {
  result <- calib.t(2450, 50)
  expect_length(result, 2)
  expect_equal(result$text$x, 2130.5, tolerance = .1)
})

### sets

test_that("pool returns average", {
  data(shroud)
  Zu <- grep("ETH", shroud$ID)
  result <- pool(shroud$y[Zu],shroud$er[Zu], talk=FALSE)
  expect_length(result, 2)
  expect_equal(result, c(676.1406, 23.7443), tolerance = 1e-5)
})

test_that("as.one returns average", {
  data(shroud)
  Zu <- grep("ETH", shroud$ID)
  result <- as.one(shroud$y[Zu],shroud$er[Zu], talk=FALSE)
  expect_equal(nrow(result), 411)
})


test_that("overlapping returns percentage", {
  y <- c(3820, 4430)
  er <- c(40, 40)
  result <- overlap(y, er, cc=1:2, talk=FALSE)
  expect_equal(result, 30.46924, tolerance=1e-5)
})

### sources

test_that("fractions returns ...", {
  Cs <- c(.02, .05, .03, .04) # carbon contents of each fraction, out of 1
  wghts <- c(5, 4, 2, .5) # weights for all fractions, e.g., in mg
  ages <- c(130, 120, 125, NA) # ages of all fractions. The unmeasured one is NA
  errors <- c(10, 12, 10, NA) # errors, unmeasured is NA
  result <- fractions(150, 20, Cs, wghts, ages, errors, talk=FALSE) # assuming a bulk age of 150 +- 20 C14 B

  expect_equal(unname(unlist(result)), c(640.474582, 436.327848), tolerance=1e-3)
})

test_that("contaminate returns correct updated result", {
  result <- contaminate(5000, 20, 1, .05, 1, talk=FALSE) # uses MC
  result <- round(result$obs,0)
  expect_equal(result[1], 4931, tolerance=2)
  expect_equal(result[2], 20, tolerance=2)
})

test_that("contaminate returns correct updated result", {
  result <- clean(5000, 20, 1, 0, 1, talk=FALSE) # no MC
  result <- round(result$obs,0)
  expect_equal(result[1], 5070.3, tolerance=2)
  expect_equal(result[2], 20, tolerance=2)
})

test_that("Monte Carlo output converges with increasing samples", {
  r1 <- clean(4000, 30, 10, 2, 1, 0.01, MC=TRUE, its=100, visualise=FALSE, talk=FALSE)
  r2 <- clean(4000, 30, 10, 2, 1, 0.01, MC=TRUE, its=5000, visualise=FALSE, talk=FALSE)

  expect_true(abs(r1$F14C["mean"] - r2$F14C["mean"]) < 0.01)
})

test_that("muck returns percentage", {
  result <- muck(600, 30, 2000, 0, 1, .01, talk=FALSE)
  expect_equal(result$perc, 67.4, tolerance=0.01)
})
