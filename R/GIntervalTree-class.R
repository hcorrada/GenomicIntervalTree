#' class definition
#' @export
setClass("GIntervalTree",
         contains="GenomicRanges",
         representation(
           seqnames="Rle",
           ranges="IRanges",
           intervalTrees="IntervalTreeList",
           rangeMap="IRangesList",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo"),
         prototype(
           seqnames=Rle(factor()),
           strand=Rle(strand()))
)

#' accessors
#' @export
setMethod("seqnames", "GIntervalTree", function(x) x@seqnames)

#' @export
setMethod("ranges", "GIntervalTree", function(x) x@ranges)

#' @export
setMethod("strand", "GIntervalTree", function(x) x@strand)

#' @export
setMethod("seqinfo", "GIntervalTree", function(x) x@seqinfo)

#' construct from GRanges object via coercion
#' @export
setAs("GRanges", "GIntervalTree",
      function(from) {
        seqnames=seqnames(from)
        rangeMap=splitRanges(seqnames)
        
        intervalTrees=IntervalTreeList(as(from, "RangesList"))
        out=new("GIntervalTree",
                seqnames=seqnames(from),
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                intervalTrees=intervalTrees,
                ranges=ranges(from),
                rangeMap=rangeMap)
        out
      })
