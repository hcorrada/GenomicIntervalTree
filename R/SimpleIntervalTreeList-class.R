#' class definition
#' 
#' @name SimpleIntervalTreeList-class
#' @family IntervalTreeList
#' 
#' @exportClass SimpleIntervalTreeList
setClass("SimpleIntervalTreeList",
         prototype=prototype(elementType="IntervalTree"),
         contains=c("IntervalTreeList", "SimpleRangesList"))

