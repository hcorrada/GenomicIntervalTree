#' findOverlaps in a PartitionedIntervalTree
#' 
#' a partition factor must be passed
#' 
#' @name findOverlaps,Ranges,PartitionedIntervalTree-method
#' @family PartitionedIntervalTree
#' 
#' @export
setMethod("findOverlaps", c("Ranges", "PartitionedIntervalTree"),
          function(query, subject,
                   maxgap=0L, minoverlap=1L,
                   type=c("any","start","within","equal"),
                   select=c("all","first","last","arbitrary"),
                   partition=NULL) {
            
            if (!isSingleNumber(maxgap) || maxgap < 0L)
              stop("'maxgap' must be a non-negative integer")
            type <- match.arg(type)
            select <- match.arg(select)
            
            if (missing(partition) || is.null(partition)) {
              stop("'partition' factor must be specified")
            }
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
            
            q_len <- length(query)
            s_len <- length(subject)
            
            s_trees <- subject@intervalTrees
            
            levelsFun <- selectMethod("levels", "Rle")
            q_levels <- levelsFun(partition)
            s_levels <- names(s_trees)            
            
            q_splitranges <- splitRanges(partition)
            
            common_levels <- intersect(q_levels, s_levels)
            results <- lapply(common_levels,
                              function(level)
                              {
                                q_idx <- q_splitranges[[level]]
                                
                                ## subjectHits are already indexed to the correspondign GRanges object
                                ## no need to remap them
                                hits <- findOverlaps(seqselect(query, q_idx), s_trees[[level]],
                                                     maxgap=maxgap,minoverlap=minoverlap,type=type)
                                
                                
                                ## remap the query hits
                                q_hits <- as.integer(q_idx)[queryHits(hits)]
                                s_hits <- subjectHits(hits)
                                not_dup <- !duplicatedIntegerPairs(q_hits,s_hits)
                                q_hits <- q_hits[not_dup]
                                s_hits <- s_hits[not_dup]
                                oo <- orderIntegerPairs(q_hits,s_hits)
                                
                                
                                new("Hits", queryHits=q_hits[oo], subjectHits=s_hits[oo],
                                    queryLength=q_len, subjectLength=s_len)
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
            
            ## hack for wierd dispatch error, it shouldn't be necessary if 
            ## code is part of GenomicRanges pacakge
            ## TODO: fix this hack 
            mergeFun=selectMethod("merge",c("Seqinfo","Seqinfo"))
            seqinfo <- mergeFun(seqinfo(query), seqinfo(subject))
            
            hits <- findOverlaps(ranges(query), ranges(subject),
                                    maxgap=maxgap,minoverlap=minoverlap,
                                    type=type,select="all",
                                    partition=seqnames(query))
             
            if (!ignore.strand) {
              q_strand <- GenomicRanges:::.strandAsSignedNumber(strand(query))
              s_strand <- GenomicRanges:::.strandAsSignedNumber(strand(subject))
              
              compatible_strand <- q_strand[queryHits(hits)] *
                s_strand[subjectHits(hits)] != -1L
              
              hits <- hits[compatible_strand]
            }
            
            q_hits <- queryHits(hits)
            s_hits <- subjectHits(hits)
            q_len <- length(query)
            s_len <- length(subject)
            
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

# #' subset a GIntervalTree object by overlaps
# #' returns a GRanges object not another GIntervalTree object
# #' 
# #' @name subsetByOverlaps,GIntervalTree,GRanges-method
# #' @family GIntervalTree
# #' @export
# #' @importClassesFrom GenomicRanges GenomicRanges
# setMethod("subsetByOverlaps", c("GIntervalTree", "GenomicRanges"),
#           function(query, subject, maxgap=0L, minoverlap=1L,
#                    type=c("any","start","end","within","equal"),
#                    ignore.strand=FALSE) {
#             hits=findOverlaps(subject, query, maxgap, minoverlap, type, select="all", ignore.strand)
#             hits=unique(subjectHits(hits))
#             as(query[hits,], "GRanges")
#           }
# )
