test_that("calBPtoBCAD returns expected BC/AD value", {
  result <- calBPtoBCAD(2000)
  expect_length(result, 1)
  expect_equal(result, -51, tolerance = 1e-5)
})

test_that("calBPtoBCAD returns expected cal BC/AD values (multiple)", {
  result <- calBPtoBCAD((-10):10)
  expect_length(result, 21)
  expect_equal(result, 1960:1940, tolerance = 1e-5)
})

test_that("calBPtob2k returns expected b2k value", {
  result <- calBPtob2k(-50)
  expect_length(result, 1)
  expect_equal(result, 0, tolerance = 1e-5)
})

test_that("calBPtoC14 returns expected C14 values", {
  result <- unlist(unname(calBPtoC14(100)))
  expect_length(result, 2)
  expect_equal(result, c(124,10), tolerance = 1e-5)
})

test_that("calBPtoF14C returns expected F14C values", {
  result <- unlist(unname(calBPtoF14C(100)))
  expect_length(result, 2)
  expect_equal(result, c(0.98468220, 0.00122656), tolerance = 1e-5)
})

test_that("calBPtopMC returns expected pMC values", {
  result <- unlist(unname(calBPtopMC(100)))
  expect_length(result, 2)
  expect_equal(result, c(98.468220, 0.122656), tolerance = 1e-5)
})

test_that("calBPtoDelta14C returns expected Delta14C values", {
  result <- unlist(unname(calBPtoDelta14C(1000)))
  expect_length(result, 2)
  expect_equal(result, c(-19.020481, 1.588828), tolerance = 1e-5)
})

###

test_that("BCADtocalBP returns expected calBP value", {
  result <- unlist(unname(BCADtocalBP(2000)))
  expect_length(result, 1)
  expect_equal(result, c(-50), tolerance = 1e-5)
})

test_that("BCADtob2k returns expected b2k value", {
  result <- unlist(unname(BCADtob2k(2000)))
  expect_length(result, 1)
  expect_equal(result, c(0), tolerance = 1e-5)
})

test_that("BCADtoC14 returns expected C14 age", {
  result <- unlist(unname(BCADtoC14(1000)))
  expect_length(result, 2)
  expect_equal(result, c(1028, 10), tolerance = 1e-5)
})

test_that("BCADtoF14C returns expected F14C value", {
  result <- unlist(unname(BCADtoF14C(1000)))
  expect_length(result, 2)
  expect_equal(result, c(0.8798779, 0.00109601), tolerance = 1e-5)
})

test_that("BCADtopMC returns expected pMC value", {
  result <- unlist(unname(BCADtopMC(1000)))
  expect_length(result, 2)
  expect_equal(result, c(87.98779, 0.109601), tolerance = 1e-5)
})

test_that("BCADtoDelta14C returns expected Delta14C value", {
  result <- unlist(unname(BCADtoDelta14C(1000)))
  expect_length(result, 2)
  expect_equal(result, c(-12.967612, 1.229487), tolerance = 1e-5)
})

###

test_that("b2ktocalBP returns expected calBP value", {
  result <- unlist(unname(b2ktocalBP(5000)))
  expect_length(result, 1)
  expect_equal(result, c(4950), tolerance = 1e-5)
})

test_that("b2ktoBCAD returns expected BC/AD value", {
  result <- b2ktoBCAD(5000)
  expect_length(result, 1)
  expect_equal(result, c(-3001), tolerance = 1e-5)
})

test_that("b2ktoC14 returns expected C14 age", {
  result <- unlist(unname(b2ktoC14(5000)))
  expect_length(result, 2)
  expect_equal(result, c(4364, 16), tolerance = 1e-5)
})

test_that("b2ktoF14C returns expected F14C value", {
  result <- unlist(unname(b2ktoF14C(5000)))
  expect_length(result, 2)
  expect_equal(result, c(0.58085213, 0.00115808), tolerance = 1e-5)
})

test_that("b2ktopMC returns expected pMC value", {
  result <- unlist(unname(b2ktopMC(5000)))
  expect_length(result, 2)
  expect_equal(result, c(58.085213, 0.115808), tolerance = 1e-5)
})

test_that("b2ktoDelta14C returns expected Delta14C value", {
  result <- unlist(unname(b2ktoDelta14C(5000)))
  expect_length(result, 2)
  expect_equal(result, c(57.077, 2.108), tolerance = 1e-3)
})

###

test_that("C14tocalBP returns expected cal BP values", {
  result <- unlist(unname(C14tocalBP(800)))
  expect_length(result, 3)
  expect_equal(result, c(690.0040, 710.0000, 717.9961), tolerance = 1e-5)
})

test_that("C14toF14C returns expected F14C value", {
  result <- C14toF14C(8000)[[1]]
  expect_type(result, "double")
  expect_length(result, 1)
  expect_equal(result, c(0.3693938), tolerance = 1e-5)
})

test_that("C14toF14C returns expected F14C values", {
  result <- unlist(unname(C14toF14C(8000, 40)))
  expect_length(result, 2)
  expect_equal(result, c(0.3693938, 0.00184397), tolerance = 1e-5)
})

test_that("C14toF14C returns expected F14C values", {
  result <- unname(C14toF14C(c(8000, 8500), c(40,35)))
  expect_length(result, 2)
  expect_equal(unlist(result), c(0.3693938, 0.3471025, 0.00184397, 0.0015163), tolerance = 1e-5)
})

test_that("F14CtoC14 returns expected C14 age", {
  result <- unname(F14CtoC14(.5, .005))
  expect_length(result, 2)
  expect_equal(unlist(result), c(5568.051301, 80.73435), tolerance = 1e-5)
})

test_that("F14CtoDelta14C returns expected Delta14C values", {
  result <- unname(F14CtoDelta14C(.5, .005, t=5500))
  expect_length(result, 2)
  expect_equal(unlist(result), c(-27.43918, 9.72561), tolerance = 1e-5)
})

test_that("Delta14CtoC14 returns expected C14 ages", {
  result <- unname(Delta14CtoC14(20, .05, t=5500))
  expect_length(result, 2)
  expect_equal(unlist(result), c(5185.477355, 0.393784), tolerance = 1e-5)
})

test_that("calBPtoF14C returns expected F14C values", {
  # Example test input
  calBP <- 2000

  # Run the function with known input
  result <- unlist(unname(calBPtoF14C(calBP)))

  # Check that the result is a numeric vector of length 2 (F14C, uncertainty)
  expect_length(result, 2)

  # Check that the values are within expected bounds
  expect_gt(result[1], 0) # F14C > 0

  # check with expected result
  expected <- c(0.772837061, 0.001348083)
  expect_equal(unlist(round(result,4)), unlist(round(expected,4)), tolerance = 1e-5)
})

