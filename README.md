A proposal for IntervalTrees on GenomicRanges objects
=================================================================

I'm writing an application where lots `findOverlap` calls are made on static `GRanges` objects. For `IRanges` we can create persistent `IntervalTree` objects that would serve the multiple overlap query use-case. There is no equivalent for `GenomicRanges` objects, so I'm proposing an implementation for this.

Please read on in the github project wiki: [`https://github.com/hcorrada/GenomicIntervalTree/wiki`](https://github.com/hcorrada/GenomicIntervalTree/wiki)

