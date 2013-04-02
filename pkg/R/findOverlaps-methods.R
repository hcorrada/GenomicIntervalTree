#' findOverlaps in an GIntervalTree
#' 
#' @name findOverlaps,GRanges,GIntervalTree-method
#' @family GIntervalTree
#' 
#' @export
#' @importClassesFrom GenomicRanges GenomicRanges
setMethod("findOverlaps", c("GenomicRanges", "GIntervalTree"),
          function(query, subject, maxgap=0L, minoverlap=1L,
                   type=c("any","start","end","within","equal"),
                   select=c("all","first","last","arbitrary"),
                   ignore.strand=FALSE) {
            if (!isSingleNumber(maxgap) || maxgap < 0L)
              stop("'maxgap' must be a non-negative integer")
            type <- match.arg(type)
            select <- match.arg(select)
            
            
            ## merge() also checks that 'query' and 'subject' are based on the
            ## same reference genome.
            
            ## hack for wierd dispatch error
            ## TODO: fix this hack
            mergeFun=selectMethod("merge",c("Seqinfo","Seqinfo"))
            seqinfo <- mergeFun(seqinfo(query), seqinfo(subject))
            
            q_len <- length(query)
            s_len <- length(subject)
            q_seqnames <- seqnames(query)
            s_seqnames <- seqnames(subject)
            
            ## another hack due to weird dispatch
            ## TODO: fix this hack
            levelsFun <- selectMethod("levels", "Rle")
            q_seqlevels <- levelsFun(q_seqnames)
            s_seqlevels <- levelsFun(s_seqnames)
            
            q_splitranges <- splitRanges(q_seqnames)
            s_splitranges <- rangeMap(subject)
            
            q_ranges <- unname(ranges(query))
            s_trees <- intervalTrees(subject)
            
            if (ignore.strand) {
              q_strand <- rep.int(1L, q_len)
              s_strand <- rep.int(1L, s_len)
            } else {
              q_strand <- GenomicRanges:::.strandAsSignedNumber(strand(query))
              s_strand <- GenomicRanges:::.strandAsSignedNumber(strand(subject))
            }
            
            common_seqlevels <- intersect(q_seqlevels, s_seqlevels)
            results <- lapply(common_seqlevels,
                              function(seqlevel)
                              {
#                                 if (isCircular(seqinfo)[seqlevel] %in% TRUE) {
#                                   circle.length <- seqlengths(seqinfo)[seqlevel]
#                                 } else {
#                                   circle.length <- NA
#                                 }
                                
                                q_idx <- q_splitranges[[seqlevel]]
                                s_idx <- s_splitranges[[seqlevel]]
#                                 hits <- .findOverlaps.circle(circle.length,
#                                                              seqselect(q_ranges, q_idx),
#                                                              seqselect(s_ranges, s_idx),
#                                                              maxgap, minoverlap, type)
                                hits <- findOverlaps(seqselect(q_ranges, q_idx), s_trees[[seqlevel]],
                                                     maxgap,minoverlap,type)
                                
                                q_hits <- queryHits(hits)
                                s_hits <- subjectHits(hits)
                                compatible_strand <- seqselect(q_strand, q_idx)[q_hits] *
                                  seqselect(s_strand, s_idx)[s_hits] != -1L
                                hits <- hits[compatible_strand]
                                remapHits(hits, query.map=as.integer(q_idx),
                                          new.queryLength=q_len,
                                          subject.map=as.integer(s_idx),
                                          new.subjectLength=s_len)
                              })
            ## Combine the results.
            q_hits <- unlist(lapply(results, queryHits))
            if (is.null(q_hits))
              q_hits <- integer(0)
            
            s_hits <- unlist(lapply(results, subjectHits))
            if (is.null(s_hits))
              s_hits <- integer(0)
            
            if (select == "arbitrary") {
              ans <- rep.int(NA_integer_, q_len)
              ans[q_hits] <- s_hits
              return(ans)
            }
            if (select == "first") {
              ans <- rep.int(NA_integer_, q_len)
              oo <- IRanges:::orderIntegerPairs(q_hits, s_hits, decreasing=TRUE)
              ans[q_hits[oo]] <- s_hits[oo]
              return(ans)
            }
            oo <- IRanges:::orderIntegerPairs(q_hits, s_hits)
            q_hits <- q_hits[oo]
            s_hits <- s_hits[oo]
            if (select == "last") {
              ans <- rep.int(NA_integer_, q_len)
              ans[q_hits] <- s_hits
              return(ans)
            }
            new2("Hits", queryHits=q_hits, subjectHits=s_hits,
                 queryLength=q_len, subjectLength=s_len,
                 check=FALSE)
          })

#' subset a GIntervalTree object by overlaps
#' returns a GRanges object not another GIntervalTree object
#' 
#' @name subsetByOverlaps,GIntervalTree,GRanges-method
#' @family GIntervalTree
#' @export
#' @importClassesFrom GenomicRanges GenomicRanges
setMethod("subsetByOverlaps", c("GIntervalTree", "GenomicRanges"),
          function(query, subject, maxgap=0L, minoverlap=1L,
                   type=c("any","start","end","within","equal"),
                   ignore.strand=FALSE) {
            hits=findOverlaps(subject, query, maxgap, minoverlap, type, select="all", ignore.strand)
            hits=unique(subjectHits(hits))
            as(query[hits,], "GRanges")
          }
)