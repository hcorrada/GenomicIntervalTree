context("GIntervalTree")

test_that("constructor works", {
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
})

test_that("constructor works 2", {
  gr <- GRanges(seqnames="chr1", ranges=IRanges(start=1:10, width=100))
  git=as(gr,"GIntervalTree")
  
  expect_equal(seqinfo(git), seqinfo(gr))
  expect_equal(seqnames(git), seqnames(gr))
  expect_equal(strand(git), strand(gr))
})

test_that("overlap works", {
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
  
  git=as(gr,"GIntervalTree")
  olaps1=findOverlaps(gr1, git)
  olaps2=findOverlaps(gr1, gr)
  expect_equal(olaps1,olaps2)
})