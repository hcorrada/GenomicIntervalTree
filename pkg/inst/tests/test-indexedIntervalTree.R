context("indexed interval tree")

test_that("indexed interval tree works", {
  x <- IRanges(start=c(1,2,3), end=c(5,2,8))
  indexes <- as.integer(10:12)
  
  z <- GenomicIntervalTree:::.indexedIntervalTree(x, indexes)
  expect_is(z, "IntervalTree")
  expect_equal(x, as(z, "IRanges"))
  expect_equal(indexes, GenomicIntervalTree:::.indexedIntervalTree_indexPositions(z))
})