![logo](https://user-images.githubusercontent.com/1739/29678184-da1f384c-88ba-11e7-9d98-df4fe3a59924.png "logo")

fast BAM/CRAM depth calculation for **WGS**, **exome**, or **targetted sequencing**.

`mosdepth` can output per-base depth about twice as fast `samtools depth`--about 25 minutes of CPU time for a 30X
genome.

it can output mean per-window depth given a window size--as would be used for CNV calling.

it can output the mean per-region given a BED file of regions.

it can create a distribution of proportion of bases covered at or above a given threshhold.

## usage

```
mosdepth

  Usage: mosdepth [options] <BAM-or-CRAM>

  -t --threads <threads>     number of CRAM|BAM decompression threads to use [default: 0]
  -c --chrom <chrom>         chromosome to restrict depth calculation.
  -Q --mapq <mapq>           mapping quality threshold [default: 0]
  -b --by <bed|window>       BED file of regions or an (integer) window-size for which to calculate depth.
  -f --fasta <fasta>         fasta file for use with CRAM files.
  -d --distribution <file>   write a cumulative distribution file (coverage, proportion).
  -h --help                  show help
```

See the section below for more info on distribution.

### exome example

To calculate the coverage in each exome capture region:
```
mosdepth --by capture.bed sample.exome.bam > sample.exome.coverage.bed
```
For a 5.5GB exome file and all 1,195,764 ensembl exons as the regions,
this completes in 1 minute 38 seconds with a single CPU.

To calculate per-base coverage:

```
mosdepth sample.exome.bam > sample.coverage.txt
```

### WGS example

For per-base whole-genome coverage:

```
mosdepth $sample.wgs.bam > $sample.txt
```

For 500-base windows (and a coverage distribution):

```
mosdepth -d $sample.dist --by 500 $sample.wgs.bam > $sample.500.bed
```

## Installation

Unless you want to install [nim](https://nim-lang.org), simply download the
[binary from the releases](https://github.com/brentp/mosdepth/releases).

If you get an error about "`libhts.so` not found", set `LD_LIBRARY_PATH`
to the directory that contains `libhts.so`. e.g.

```LD_LIBRARY_PATH=~/src/htslib/ mosdepth -h```

## distribution output

This is **useful for QC**.

The `--distribution` option outputs a cumulative distribution indicating then
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

Using something like that expanded for use on multiple samples, we can
plot the distribution from the entire genome. Below we show this for samples
with ~60X coverage:

![WGS Example](https://user-images.githubusercontent.com/1739/29646192-2a2a6126-883f-11e7-91ab-049295eb3531.png "WGS Example")

We can also run `mosdepth` on just the Y chromosome (--chrom Y) to verify that males and females
track separately. Below, we that see female samples cluster along the axes while male samples have
close to 30X coverage for almost 40% of the genome.

![Y Example](https://user-images.githubusercontent.com/1739/29646191-2a246564-883f-11e7-951a-aa68d7a1a6ed.png "Y Example")

See [this blog post](http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html) for
more details.

## how it works

As it encounters each chromosome, `mosdepth` creates an array the length of the chromosome.
For every start it encounters, it increments the value in that position of the array. For every
stop, it decrements that position. From this, the depth at a particular position is the
cumulative sum of all array positions preceding it (a similar algorithm is used in BEDTools
where starts and stops are tracked separately). `mosdepth` **avoids double-counting
overlapping mate-pairs** and it **tracks every aligned part of every read using the CIGAR
operations**. Because of this data structure, the the coverage `distribution` calculation
can be done without a noticeable increase in run-time. The image below conveys the concept:

![alg](https://user-images.githubusercontent.com/1739/29647913-d79ab028-8848-11e7-86cf-60d4b087bc3b.png "algorithm")

This array accounting is very fast. There are no extra allocations or objects to track and
it is also conceptually simple. For these reasons, it is faster than `samtools depth` which
works by using the [pileup](http://samtools.sourceforge.net/pileup.shtml) machinery that
tracks each read, each base. 

The `mosdepth` method has some limitations. Because a large array is allocated and it is
required (in general) to take the cumulative sum of all preceding positions to know the depth
at any position, it is slower for small, 1-time regional queries. It is, however fast for
window-based or BED-based regions, because it first calculates the full chromosome coverage
and then reports the coverage for each region in that chromosome. Another downside is it uses
more memory than samtools. The amount of memory is approximately equal to 32-bits * longest chrom
length, so for the 249MB chromosome 1, it will require 500MB of memory.

`mosdepth` is written in [nim](https://nim-lang.org/) and it uses our [htslib](https://github.com/samtools/htslib)
via our nim wrapper [hts-nim](https://github.com/brentp/hts-nim/)

## output

When the `--by` argument is *not* specified, the output of `mosdepth`. The output looks like:

```
chr1	216	0
chr1	232	1
chr1	236	2
chr1	255	1
chr1	9991	0
chr1	9992	1
chr1	9993	2
chr1	9994	6
chr1	9995	5
chr1	9996	6
```

Each line indicates the end of the previous coverage level and the start of the next.

So the first line indicates that:

+ the values on chr1 from 0 to 216 (in 0-based, half-open (BED) coordinates) have a depth of 0. 
+ the values on chr1 from 216 to 232 have a depth of 1 
+ 232..236 == 2
+ 236..255 == 1
+ 255..9991 == 0
+ and so on ...

## speed and comparison

coming soon.

format: BAM/CRAM
tool: samtools, sambamba, bedtools, mosdepth
threads: 0..5
mode: base | window

----------------------------------------------------------------------------
 format |    tool    | threads  | mode | relative speed | run-time | memory |


o
# plot of per-base samtools vs per-base mosdepth
