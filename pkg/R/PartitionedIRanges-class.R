#' class defintion
#' 
#' @name PartitionedIRanges-class
#' @family PartitionedIRanges
#' 
#' @exportClass PartitionedIRanges
setClass("PartitionedIRanges",
         contains="IRanges",
         representation(partition="Rle"))

#' constructor
#' 
#' @rdname PartitionedIRanges-class
#' @family PartitionedIRanges
#' @export
PartitionedIRanges <- function(ranges, partition) {
  
  out <- new2("PartitionedIRanges", start=start(ranges), width=width(ranges), NAMES=names(ranges), partition=partition, check=FALSE)
}

#' partition accessor generic
#'
#' @rdname PartitionedIRanges-class
#' @family PartitionedIRanges 
#' @export
setGeneric("partition", function(x) standardGeneric("partition"))

#' partition accessor method
#' 
#' @rdname PartitionedIRanges-class
#' @family PartitionedIRanges
#' @export
setMethod("partition", "PartitionedIRanges", function(x) x@partition)

#' subset method
#' 
#' @rdname PartitionedIRanges-class
#' @family PartitionedIRanges
#' @export
setMethod("[", "PartitionedIRanges", 
          function(x, i, j, ...) {
            out <- callNextMethod()
            slot(out, "partition", check=FALSE) <- x@partition[i]
            out
          })
