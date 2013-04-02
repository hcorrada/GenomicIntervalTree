A proposal for persistent IntervalTrees on GenomicRanges objects
=================================================================

I'm writing an application where many `findOverlap` calls are made on relatively static `GRanges` objects. For `IRanges` we can create persistent `IntervalTree` objects that would serve the multiple overlap query use-case. There is no equivalent for `GenomicRanges` objects, so I'm proposing this.

Usage
------

This is the intended usage:

```{r}
gr <-
    GRanges(seqnames =
              Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges =
              IRanges(1:10, width = 10:1, names = head(letters,10)),
            strand =
              Rle(strand(c("-", "+", "*", "+", "-")),
                  c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10))
  gr1 <-
    GRanges(seqnames = "chr2", ranges = IRanges(4:3, 6),
            strand = "+", score = 5:4, GC = 0.45)
  
  git <- GIntervalTree(gr)
  
  # this uses the interval trees for each seqlevel
  olaps1 <- findOverlaps(gr1, git)
```
Classes
--------

### IntervalTreeList

This class defines a list of `IntervalTree` objects using the `List` interface defined in `IRanges`. We split the ranges in `GenomicRanges` objects by `seqlevel` and create an IntervalTree for each range. We store those trees in an `IntervalTreeList` object.

I also defined `SimpleIntervalTreeList` as an implementation of this VIRTUAL class.

### GIntervalTree

This is the main class defined. It extends `GenomicRanges` so its interface is supported.
Given a `GenomicRanges` object, a list of IntervalTrees are constructed for the ranges for each seqlevel and stored in the `intervalTrees` slot. A vector mapping indices of the original ranges to the corresponding IntervalTree is also stored.

Methods
--------

### coerce

A number of `coerce` methods are defined. Main ones are `coerce,RangesList,IntervalTreeList-method` and `coerce,GenomicRanges,GIntervalTree-method`. 

### findOverlaps

The `findOverlaps,GenomicRanges,GIntervalTree-method` is defined here. For each of the common seqlevels in the  query and the subject, the corresponding IntervalTree of the subject is used. All other derived methods, e.g., `countOverlaps` and `subsetByOverlaps` where the subject is a `GIntervalOverlaps` is the subject call this method.

### subsetByOverlaps
**THIS IS WIERD!** The `subsetByOverlaps,GIntervalTree,GenomicRanges-method` is also implemented here. It returns a `GRanges` object (instead of an `GIntervalTree` object). The expectation is that objects `sgr1` and `sgr2` in the following example are equal.

```{r}
sgr1 <- subsetByOverlaps(query, subject)
sgr2 <- subsetByOverlaps(GIntervalTree(query), subject)
```



