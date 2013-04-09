#' union class
#' @exportClass IRangesOrRangesList
#' @import GenomicRanges
setClassUnion("IRangesOrRangesList", c("IRanges", "RangesList"))

#' GIntervalTree class
#' 
#' Defines persistent interval trees for GRanges objects.
#' 
#' 
#' @name GIntervalTree-class
#' @family GIntervalTree
#' 
#' @exportClass GIntervalTree
#' @import GenomicRanges
#' @import BiocGenerics
setClass("GIntervalTree",
         contains="GenomicRanges",
         representation(
           seqnames="Rle",
           ranges="IRangesOrRangesList",
           #intervalTrees="IntervalTreeList",
           #rangeMap="IRangesList",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo"),
         prototype(
           seqnames=Rle(factor()),
           strand=Rle(strand()))
)

#' seqnames accessor 
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges seqnames
setMethod("seqnames", "GIntervalTree", function(x) x@seqnames)

#' ranges accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges ranges
setMethod("ranges", "GIntervalTree", function(x) x@ranges)

#' strand accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges strand
setMethod("strand", "GIntervalTree", function(x) x@strand)

#' seqinfo accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges seqinfo
setMethod("seqinfo", "GIntervalTree", function(x) x@seqinfo)

#' intervalTrees accessor generic
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
setGeneric("intervalTrees", function(x) standardGeneric("intervalTrees"))

#' intervalTrees accessor method
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
setMethod("intervalTrees", "GIntervalTree", function(x) x@intervalTrees)

#' rangesMap accessor generic
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
setGeneric("rangeMap", function(x) standardGeneric("rangeMap"))

#' rangesMap accessor method
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
setMethod("rangeMap", "GIntervalTree", function(x) x@rangeMap)


#' construct from GRanges object via coercion
#' 
#' @name as
#' @family GIntervalTree
#' @importClassesFrom GenomicRanges GRanges
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
                #intervalTrees=intervalTrees,
                #ranges=ranges(from),
                ranges=intervalTrees)
                #rangeMap=rangeMap)
        out
      }
)

#' constructor function using GRanges object
#' 
#' @family GIntervalTree
#' @export
GIntervalTree <- function(x) {
  as(x, "GIntervalTree")
}

#' coercion from GIntervalTree to GRanges object
#' 
#' @family GIntervalTree
#' @name as
#' @importClassesFrom GenomicRanges GRanges
setAs("GIntervalTree", "GRanges",
      function(from) {
        out=new("GRanges",
                seqnames=seqnames(from),
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                ranges=ranges(from))
      }
)
