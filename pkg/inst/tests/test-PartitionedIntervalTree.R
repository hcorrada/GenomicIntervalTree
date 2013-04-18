context("PartitionedIntervalTree")

test_that("constructor works", {
  x1 <- IRanges(start=c(1,2,3), end=c(5,2,8))
  x2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  x <- c(x1,x2)
  partition <- Rle(factor(rep(c("one","two"),c(3,4))))
  
  y <- PartitionedIntervalTree(x, partition)
  
  expect_is(y,"PartitionedIntervalTree")
  expect_equal(x, as(y, "IRanges"))
})

test_that("subset works", {
  x1 <- IRanges(start=c(1,2,3), end=c(5,2,8))
  x2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  x <- c(x1,x2)
  partition <- Rle(factor(rep(c("one","two"),c(3,4))))
  
  y <- PartitionedIntervalTree(x, partition)
  y2 <- y[2:4]
  
  expect_is(y2,"PartitionedIntervalTree")
  expect_equal(x[2:4], as(y2, "IRanges"))
  expect_equal(GenomicIntervalTree:::.getPartition(y2), partition[2:4])
})


test_that("findOverlaps works", {
  subject <- IRanges(start=c(1,2,3,15,45,20,1), end=c(5,2,8,15,100,80,5))
  s_partition <-Rle(factor(rep(c("one","two"),c(3,4))))
  tsubject <- PartitionedIntervalTree(subject, s_partition)
  
  query <- IRanges(start=c(1,3,40), end=c(10,20,120))
  q_partition <- Rle(factor(c("one","one","two")))

  olaps2 <- findOverlaps(query, tsubject, partition=q_partition)
  
  for (lvl in levels(q_partition)) {
    s_idx <- which(as.logical(s_partition == lvl))
    q_idx <- which(as.logical(q_partition == lvl))
    
    olaps1 <- findOverlaps(query[q_idx],subject[s_idx])
    expect_equal(q_idx[queryHits(olaps1)],queryHits(olaps2)[queryHits(olaps2) %in% q_idx])
    expect_equal(s_idx[subjectHits(olaps1)], subjectHits(olaps2)[subjectHits(olaps2) %in% s_idx])
  }
})