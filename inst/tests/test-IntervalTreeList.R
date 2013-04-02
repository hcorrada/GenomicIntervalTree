context("IntervalTreeList")

test_that("coerce works", {
  x <- IRangesList(start=list(c(1,2,3), c(15,45,20,1)),
                   end=list(c(5,2,8), c(15,100,80,5)))
  y <- as(x,"IntervalTreeList")
  
  expect_is(y,"SimpleIntervalTreeList")
  expect_equal(x[[1]], as(y[[1]], "IRanges"))
  expect_equal(x[[2]], as(y[[2]], "IRanges"))
})

test_that("constructor function works", {
  x <- IRangesList(start=list(c(1,2,3), c(15,45,20,1)),
                   end=list(c(5,2,8), c(15,100,80,5)))
  y <- IntervalTreeList(x)
  
  expect_is(y,"SimpleIntervalTreeList")
  expect_equal(x[[1]],as(y[[1]],"IRanges"))
  expect_equal(x[[2]], as(y[[2]], "IRanges"))
  
})

test_that("findOverlaps works", {
  x <- IRangesList(start=list(c(1,2,3), c(15,45,20,1)),
                   end=list(c(5,2,8), c(15,100,80,5)))
  tx <- IntervalTreeList(x)
  
  y <- IRangesList(start=list(c(1,3), c(40)), end=list(c(10,20),c(120)))
  
  olaps1 <- findOverlaps(y, x)
  olaps2 <- findOverlaps(y, tx)
  
  expect_equal(olaps1,olaps2)
})