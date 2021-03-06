context("GIntervalTree")

test_that("coercion on fully specified GRanges works", {
  seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
  gr <-
    GRanges(seqnames =
              Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges = IRanges(
              1:10, width = 10:1, names = head(letters,10)),
            strand = Rle(
              strand(c("-", "+", "*", "+", "-")),
              c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10),
            seqinfo=seqinfo)
  git=as(gr,"GIntervalTree")
  
  expect_equal(seqinfo(git), seqinfo(gr))
  expect_equal(seqnames(git), seqnames(gr))
  expect_equal(strand(git), strand(gr))
  
  git_rng <- as(ranges(git), "IRanges")
  expect_equal(start(git_rng), start(gr))
  expect_equal(width(git_rng), width(gr))
  expect_equal(GenomicIntervalTree:::.getPartition(ranges(git)), seqnames(gr))
})

test_that("constructor on fully specified GRanges works", {
  seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
  gr <-
    GRanges(seqnames =
              Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges = IRanges(
              1:10, width = 10:1, names = head(letters,10)),
            strand = Rle(
              strand(c("-", "+", "*", "+", "-")),
              c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10),
            seqinfo=seqinfo)
  git <- GIntervalTree(gr)
  
  expect_equal(seqinfo(git), seqinfo(gr))
  expect_equal(seqnames(git), seqnames(gr))
  expect_equal(strand(git), strand(gr))
  
  git_rng <- as(ranges(git), "IRanges")
  expect_equal(start(git_rng), start(gr))
  expect_equal(width(git_rng), width(gr))
  expect_equal(GenomicIntervalTree:::.getPartition(ranges(git)), seqnames(gr))
  
})

test_that("coercion on partially defined GRanges works", {
  gr <- GRanges(seqnames=rep(c("chr1","chr2"),len=10), ranges=IRanges(start=1:10, width=100))
  git=as(gr,"GIntervalTree")
  
  expect_equal(seqinfo(git), seqinfo(gr))
  expect_equal(seqnames(git), seqnames(gr))
  expect_equal(strand(git), strand(gr))
  
  git_rng <- as(ranges(git), "IRanges")
  expect_equal(start(git_rng), start(gr))
  expect_equal(width(git_rng), width(gr))
  expect_equal(GenomicIntervalTree:::.getPartition(ranges(git)), seqnames(gr))
})

test_that("conversion to GRanges is equal (except names)", {
  gr <- GRanges(seqnames=rep(c("chr1","chr2"),len=10), ranges=IRanges(start=1:10, width=100))
  gr2 <- as(as(gr,"GIntervalTree"), "GRanges")
  
  expect_equal(gr,gr2)
})

test_that("subsetting works", {
  gr <- GRanges(seqnames=rep(c("chr1","chr2"),len=10), ranges=IRanges(start=1:10, width=100))
  git=as(gr,"GIntervalTree")
  
  gr <- gr[4:7,]
  git <- git[4:7,]
  
  expect_equal(seqinfo(git), seqinfo(gr))
  expect_equal(seqnames(git), seqnames(gr))
  expect_equal(strand(git), strand(gr))
  
  git_rng <- as(ranges(git), "IRanges")
  expect_equal(start(git_rng), start(gr))
  expect_equal(width(git_rng), width(gr))
  expect_equal(GenomicIntervalTree:::.getPartition(ranges(git)), seqnames(gr))
  
})

test_that("findOverlaps works", {
  gr <-
    GRanges(seqnames =
              Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges =
              IRanges(1:10, width = 10:1, names = head(letters,10)),
            strand =
              Rle(strand(c("-", "+", "*", "+", "-")),
                  c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10))
  gr1 <-
    GRanges(seqnames = "chr2", ranges = IRanges(4:3, 6),
            strand = "+", score = 5:4, GC = 0.45)
  
  git <- GIntervalTree(gr)
  
  olaps2=findOverlaps(gr1, gr)
  olaps1=findOverlaps(gr1, git)
  expect_equal(olaps1,olaps2)
})

# test_that("subsetByOverlaps,GenomicRanges,GIntervalTree-method works", {
#   gr <-
#     GRanges(seqnames =
#               Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#             ranges =
#               IRanges(1:10, width = 10:1, names = head(letters,10)),
#             strand =
#               Rle(strand(c("-", "+", "*", "+", "-")),
#                   c(1, 2, 2, 3, 2)),
#             score = 1:10,
#             GC = seq(1, 0, length=10))
#   gr1 <-
#     GRanges(seqnames = "chr2", ranges = IRanges(4:3, 6),
#             strand = "+", score = 5:4, GC = 0.45)
#   
#   git=as(gr,"GIntervalTree")
#   
#   olaps2=subsetByOverlaps(gr1, gr)
#   olaps1=subsetByOverlaps(gr1, git)
#   expect_equal(olaps1,olaps2)
# })
# 
# test_that("subsetByOverlaps,GIntervalTree,GenomicRanges-method works", {
#   gr <-
#     GRanges(seqnames =
#               Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#             ranges =
#               IRanges(1:10, width = 10:1, names = head(letters,10)),
#             strand =
#               Rle(strand(c("-", "+", "*", "+", "-")),
#                   c(1, 2, 2, 3, 2)),
#             score = 1:10,
#             GC = seq(1, 0, length=10))
#   gr1 <-
#     GRanges(seqnames = "chr2", ranges = IRanges(4:3, 6),
#             strand = "+", score = 5:4, GC = 0.45)
#   
#   git=as(gr,"GIntervalTree")
#   
#   olaps2=subsetByOverlaps(gr, gr1)
#   olaps1=subsetByOverlaps(git, gr1)
#   expect_equal(olaps1,olaps2)
# })