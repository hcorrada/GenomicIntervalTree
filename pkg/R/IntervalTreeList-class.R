#' class definition
#' @name IntervalTreeList-class
#' @family IntervalTreeList
#' 
#' @exportClass IntervalTreeList
setClass("IntervalTreeList",
         representation("VIRTUAL"),
         prototype=prototype(elementType="IntervalTree"),
         contains="RangesList")

#' constructor
#' @export
IntervalTreeList <- function(rangesList) {
  as(rangesList, "IntervalTreeList")
}

#' construct from ranges list by coercion
#' @name as
#' @family IntervalTreeList
#' @importClassesFrom IRanges RangesList
setAs("RangesList", "IntervalTreeList",
      function(from) {
        #IRanges:::validObject(from)
        listData=vector("list", length(from))
        for (i in seq_along(from))
          listData[[i]] <- as(from[[i]], "IntervalTree")
        names(listData) <- names(from)
        IRanges:::newList("SimpleIntervalTreeList", listData)
      }
)

