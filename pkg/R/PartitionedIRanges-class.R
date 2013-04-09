#' class defintion
#' @exportClass PartitionedIRanges
setClass("PartitionedIRanges",
         contains="IRanges",
         representation=(partition="Rle"))

#' constructor
#' @export
PartitionedIRanges <- function(ranges, partition) {
  if (class(ranges) != "IRanges")
    stop("'ranges' must be of class 'IRanges'")
  
  if (is(partition, "Rle")) {
    if (!is.factor(runValue(partition))) {
      stop("'partition' must be a 'factor' Rle or 'factor'")
    }
  } else {
    if (!is.factor(partition)) {
      stop("'partition' must be a 'factor' Rle or 'factor'")
    }
    partition <- Rle(partition)
  }
  
  out <- new2(start=start(ranges), end=end(ranges), NAMES=names(ranges), partition=partition, check=FALSE)
}
