A proposal for IntervalTrees on GenomicRanges objects
=================================================================

I'm writing an application where lots `findOverlap` calls are made on static `GRanges` objects. For `IRanges` we can create persistent `IntervalTree` objects that would serve the multiple overlap query use-case. There is no equivalent for `GenomicRanges` objects, so I'm proposing this.

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
  olaps <- findOverlaps(gr1, git)
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

A lot of the code for the `findOverlaps,GenomicRanges,GIntervalTree-method` method is duplicated from the `findOverlaps,GenomicRanges,GenomicRanges-method`. Some of this code is there to deal with splitting up the ranges in the `GenomicRanges` object. It might make more sense for these functions to call a `findOverlaps,RangesList,RangesList-method`?

### subsetByOverlaps
**THIS IS WEIRD!** The `subsetByOverlaps,GIntervalTree,GenomicRanges-method` is also implemented here. It returns a `GRanges` object (instead of an `GIntervalTree` object). The expectation is that objects `sgr1` and `sgr2` in the following example are equal.

```{r}
sgr1 <- subsetByOverlaps(query, subject)
sgr2 <- subsetByOverlaps(GIntervalTree(query), subject)
```

There is no equivalent method for `IntervalTree` objects, e.g., `subsetByOverlaps,IntervalTree,Ranges-method`. Perhaps it should be defined in a similar way. Otherwise, a new generic might be worth defining that has this functionality.
 
One last thing
---------------

The implementation is not complete yet (it just implements what I need for my application). 