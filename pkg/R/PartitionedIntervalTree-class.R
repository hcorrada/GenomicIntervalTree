#' definition
#' @exportClass PartitionedIntervalTree
#' @import IRanges
setClass("PartitionedIntervalTree",
         contains="Ranges",
         representation(intervalTrees="IntervalTreeList",
                        length="integer"))

.indexedIntervalTree <- function(ranges, indexes) {
  ptr <- .Call2("IntegerIndexedIntervalTree_new", ranges, indexes, PACKAGE="IRanges")
  new2("IntervalTree", ptr=ptr, mode="integer", check=FALSE)
}

.indexedIntervalTree_indexPositions <- function(tree) {
  .Call2("IntegerIndexedIntervalTree_indexPositions", tree@ptr, PACKAGE="IRanges")
}

#' constructor
#' @family PartitionedIntervalTree
#' @export
PartitionedIntervalTree <- function(ranges, partition) {
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
  
  indexes <- splitRanges(partition)
  intervalTrees <- vector("list", length(indexes))
  names(intervalTrees) <- names(indexes)
  for (i in seq_along(indexes)) {
    intervalTrees[[i]] <- .indexedIntervalTree(seqselect(ranges, indexes[[i]]), as.integer(indexes[[i]]))
  }
  intervalTrees <- IRanges:::newList("SimpleIntervalTreeList", intervalTrees)
  new2("PartitionedIntervalTree", intervalTrees=intervalTrees, length=length(ranges), check=FALSE)
}

#' coerce back to IRanges
#' @name as
setAs("PartitionedIntervalTree", "IRanges",
      function(from) {
        start <- integer(length(from))
        width <- integer(length(from))
        
        for (i in seq_along(from@intervalTrees)) {
          curTree <- from@intervalTrees[[i]]
          cur_ranges <- as(curTree, "IRanges")
          cur_indexes <- .indexedIntervalTree_indexPositions(curTree)
          
          start[cur_indexes] <- start(cur_ranges)
          width[cur_indexes] <- width(cur_ranges)
        }
        new2("IRanges", start=start, width=width, check=FALSE)
      })

#' length accessor
#' @family PartitionedIntervalTree
#' @export
setMethod("length", "PartitionedIntervalTree", function(x) x@length)

