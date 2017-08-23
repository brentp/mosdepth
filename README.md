# mosdepth

fast BAM/CRAM depth calculation

`mosdepth` can output per-base depth about twice as fast `samtools depth`.

it can output per-window depth given a window size or a BED file of regions.

and it can create a distribution of proportion of bases covered at or above
a given threshhold.

## usage

```
mosdepth

  Usage: mosdepth [options] <BAM-or-CRAM>

  -t --threads <threads>     number of CRAM|BAM decompression threads to use [default: 0]
  -c --chrom <chrom>         chromosome to restrict depth calculation.
  -Q --mapq <mapq>           mapping quality threshold [default: 0]
  -b --by <bed|window>       BED file of regions or an (integer) window-size for which to calculate depth.
  -f --fasta <fasta>         fasta file for use with CRAM files.
  -d --distribution <file>   write a cumulative distribution file (coverage, n).
  -h --help                  show help
```

See the section below for more info on coverage.

### exome example

To calculate the coverage in each exome capture region:
```
mosdepth -b capture.bed sample.exome.bam > sample.exome.coverage.bed
```
For a 5.5GB exome file and all 1,195,764 ensembl exons as the regions,
this completes in 1 minute 38 seconds with a single CPU.


### WGS example

For per-base whole-genome coverage:

```
mosdepth sample.wgs.bam > sample.txt
```

For 500-base windows:

```
mosdepth --by 500 sample.wgs.bam > sample.500.bed
```

## distribution

`distribution` gives the proportion of bases covered at a given threshhold.

the `--distribution` option, outputs a cumulative distribution indicating then
proportion of genome bases (or the proportion of the `--by`) that were covered
for at least a given coverage value.

This will write a file to the requested path with values indicating the coverage
threshold and the proportion of bases covered at that threshold.

This could be plotted in python with e.g.:

```
from matplotlib import pyplot as plt
xs, ys = zip(`*`(map(float, x.split()) for x in open('dist.txt')))
plt.plot(xs, ys)
```

see [this blog post](http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html) for
more details.

## how it works

For each chromosome, `mosdepth` creates an array the length of the chromosome. For every
start it encounters, integer increments value in that position of the array. For every
stop, it decrements that position. From this, the depth at a particular position is the
cumulative sum of all array positions preceding it (a similar algorithm is used in BEDTools
where starts and stops are tracked separately). `mosdepth` avoids double-counting
overlapping mate-pairs and it tracks every aligned part of every read using the CIGAR
operations. Because of this data structure, the the coverage `distribution` calculation
is also very fast.

This array accounting is very fast. There are no extra allocations or objects to track and
it is also conceptually simple. For these regions, it is faster than `samtools depth` which
works by using the [pileup](http://samtools.sourceforge.net/pileup.shtml) machinery that
tracks each read, each base. 

This method has some limitations. Because a large array is allocated and it is required (in general)
to take the cumulative sum of all preceding positions to know the depth at any position,
it is slower for small, 1-time regional queries. It is, however fast for window-based
or BED-based regions, because it first calculates the full chromosome coverage and then
reports the coverage for each region in that chromosome.

## output

When no `--by` argument is specified, the output of `mosdepth`. The output looks like:
```
```

Each line indicates the end of the previous coverage level and the start of the next.

So the first line indicates that:

+ the values on chr1 from 0 to 216 (in 0-based, half-open (BED) coordinates) have a depth of 0. 
+ the values on chr1 from 216 to 232 have a depth of 1 
+ and so on ...

