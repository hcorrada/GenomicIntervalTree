#' definition
#' @exportClass PartitionedIntervalTree
#' @import IRanges
setClass("PartitionedIntervalTree",
         representation=representation("intervalTrees"="list"),
         contains=c("IntervalTree","IRanges"))

#' constructor
#' @export
PartitionedIntervalTree(from, partition) {
  if (class(from) != "IRanges")
    stop("argument 'from' must be of class IRanges")
  
  IRanges:::validObject(from)
  
  indexes <- splitRanges(partition)
  intervalTrees <- vector("list", length(indexes))
  names(intervalTrees) <- names(indexes)
  for (i in seq_along(indexes)) {
    ptr <- .Call2("IntegerIndexedIntervalTree_new", from[indexes[[i]]], as.integer(indexes[[i]]), PACKAGE="GenomicIntervalTree")
    intervalTrees[[i]] <- new2("IntervalTree", ptr=ptr, mode="integer", check=FALSE)
  }
  
  out <- new2("PartitionedIntervalTree", intervalTrees=intervalTrees, check=FALSE)
}
