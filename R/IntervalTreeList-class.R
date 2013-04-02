#' class definition
#' @export
setClass("IntervalTreeList",
         representation("VIRTUAL"),
         prototype=prototype(elementType="IntervalTree"),
         contains="RangesList")

#' constructor
#' @export
IntervalTreeList=function(rangesList) {
  as(rangesList, "IntervalTreeList")
}

#' construct from ranges list by coercion
#' @export
setAs("RangesList", "IntervalTreeList",
      function(from) {
        #IRanges:::validObject(from)
        listData=vector("list", length(from))
        for (i in seq_along(from))
          listData[[i]] = as(x[[i]], "IntervalTree")
        IRanges:::newList("SimpleIntervalTreeList", listData)
      })
