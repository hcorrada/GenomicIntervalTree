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
           ranges="PartitionedIntervalTree",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo"),
         prototype(
           strand=Rle(strand()))
)

.valid.GIntervalTree.length <- function(x) {
  n <- length(ranges(x))
  if ((length(strand(x)) != n)
      || (nrow(mcols(x)) != n))
    return("slot lengths are not all equal")
  NULL
}

.valid.GIntervalTree.ranges <- function(x) {
  if (class(ranges(x)) != "PartitionedIntervalTree")
    return("'ranges(x)' must be a PartitionedIntervalTree instance")
  NULL
}

.valid.GIntervalTree <- function(x) {
  c(.valid.GIntervalTree.length(x),
    .valid.GIntervalTree.ranges(x),
    .valid.GenomicRanges.strand(x),
    .valid.GenomicRanges.mcols(x),
    valid.GenomicRanges.seqinfo(x))
}

setValidity2("GIntervalTree", .valid.GIntervalTree)

#' seqnames accessor 
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges seqnames
setMethod("seqnames", "GIntervalTree", function(x) .getPartition(ranges(x)))

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

#' construct from GRanges object via coercion
#' 
#' @name as
#' @family GIntervalTree
#' @importClassesFrom GenomicRanges GRanges
setAs("GRanges", "GIntervalTree",
      function(from) {
        out=new2("GIntervalTree",
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                ranges=PartitionedIntervalTree(ranges(from), seqnames(from)),
                check=FALSE)
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
                seqnames=.getPartition(ranges(from)),
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                ranges=as(ranges(from), "IRanges"))
      }
)

#' subsetting
#' 
#' @family GIntervalTree
#' @rdname GIntervalTree-class
#' @export
setMethod("[", "GIntervalTree",
          function(x, i, j, ...) {
            gr <- as(x, "GRanges")[i]
            as(gr, "GIntervalTree")
          })