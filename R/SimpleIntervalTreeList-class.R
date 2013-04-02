#' class definition
#' @export
setClass("SimpleIntervalTreeList",
         prototype=prototype(elementType="IntervalTree"),
         contains=c("IntervalTreeList", "SimpleRangesList"))

