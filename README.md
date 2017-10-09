fast BAM/CRAM depth calculation for **WGS**, **exome**, or **targeted sequencing**.

![logo](https://user-images.githubusercontent.com/1739/29678184-da1f384c-88ba-11e7-9d98-df4fe3a59924.png "logo")

[![Build Status](https://travis-ci.org/brentp/mosdepth.svg?branch=master)](https://travis-ci.org/brentp/mosdepth)

`mosdepth` can output:

+ per-base depth about 2x as fast `samtools depth`--about 25 minutes of CPU time for a 30X genome.
+ mean per-window depth given a window size--as would be used for CNV calling.
+ the mean per-region given a BED file of regions.
+ a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide.

when appropriate, the output files are bgzipped and indexed for ease of use.

## usage

```
mosdepth

**NOTE**: the usage below is for the development (coming) version of mosdepth which differs from the current release.

  Usage: mosdepth [options] <prefix> <BAM-or-CRAM>

Arguments:

  <prefix>       outputs: `{prefix}.mosdepth.dist.txt`
                          `{prefix}.per-base.bed.gz` (unless -n/--no-per-base is specified)
                          `{prefix}.regions.bed.gz` (if --by is specified)
                          `{prefix}.quantized.bed.gz` (if --quantize is specified)

  <BAM-or-CRAM>  the alignment file for which to calculate depth.

Common Options:
  
  -t --threads <threads>     number of BAM decompression threads [default: 0]
  -c --chrom <chrom>         chromosome to restrict depth calculation.
  -b --by <bed|window>       optional BED file or (integer) window-sizes.
  -n --no-per-base           dont output per-base depth (skipping this output will speed execution).
  -f --fasta <fasta>         fasta file for use with CRAM files.

Advanced options:

  -F --flag <FLAG>           exclude reads with any of the bits in FLAG set [default: 1796]
  -q --quantize <segments>    write quantized output see docs for description.
  -Q --mapq <mapq>           mapping quality threshold [default: 0]
  -h --help                  show help
```

See the section below for more info on distribution.

If `--by` is a BED file with 4 or more columns, it is assumed the the 4th column is the name.
That name will be propagated to the `mosdepth` output in the 4th column with the depth in the 5th column.
If you don't want this behavior, simply send a bed file with 3 columns.

### exome example

To calculate the coverage in each exome capture region:
```
mosdepth --by capture.bed sample-output sample.exome.bam
```
For a 5.5GB exome BAM and all 1,195,764 ensembl exons as the regions,
this completes in 1 minute 38 seconds with a single CPU.

Per-base output will go to `sample-output.per-base.bed.gz`,
the mean for each region will go to `sample-output.regions.bed.gz`;
each of those will be written along with a CSI index that can be
used for tabix queries.
The distribution of depths will go to `sample-output.mosdepth.dist.txt`

### WGS example

For 500-base windows

```
mosdepth -n --by 500 sample.wgs $sample.wgs.bam
```

`-n` means don't output per-base data, this will make `mosdepth`
a bit faster as there is some cost to outputting that much text.

### Distribution only

To get only the distribution value, without the depth file or the per-base and using 3 threads:

```
mosdepth -n -t 3 $sample $bam
```

Output will go to `$sample.mosdepth.dist.txt`

## Installation

Unless you want to install [nim](https://nim-lang.org), simply download the
[binary from the releases](https://github.com/brentp/mosdepth/releases).

`mosdepth` uses requires htslib version 1.4 or later. If you get an error 
about "`libhts.so` not found", set `LD_LIBRARY_PATH` to the directory that 
contains `libhts.so`. e.g.

`LD_LIBRARY_PATH=~/src/htslib/ mosdepth -h`

If you get the error `could not import: hts_check_EOF` you may need to 
install a more recent version of htslib.

If you do want to install from source, see the [travis.yml](https://github.com/brentp/mosdepth/blob/master/.travis.yml)
and the [install.sh](https://github.com/brentp/mosdepth/blob/master/scripts/install.sh).

## distribution output

This is **useful for QC**.

The `$prefix.mosdepth.dist.txt` file contains, a cumulative distribution indicating the
proportion of bases (or the proportion of the `--by`) that were covered
for at least a given coverage value. It does this for each chromosom, and for then
whole genome.

Each row will indicate:
 + chromosome (or "genome")
 + coverage level
 + proportion of bases covered at that level

The last value in each chromosome will be coverage level of 0 aligned with
1.0 bases covered at that level.

A python plotting script is provided in `scripts/plot-dist.py` that will make 
plots like below. Use is `python scripts/plot-dist.py *.dist` and the output
is `dist.html` with a plot for the full set along with one for each chromosome.

Using something like that, we can plot the distribution from the entire genome.
Below we show this for samples with ~60X coverage:

![WGS Example](https://user-images.githubusercontent.com/1739/29646192-2a2a6126-883f-11e7-91ab-049295eb3531.png "WGS Example")

We can also view the Y chromosome to verify that males and females
track separately. Below, we that see female samples cluster along the axes while male samples have
close to 30X coverage for almost 40% of the genome.

![Y Example](https://user-images.githubusercontent.com/1739/29646191-2a246564-883f-11e7-951a-aa68d7a1a6ed.png "Y Example")

See [this blog post](http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html) for
more details.

## quantize

quantize allows splitting coverage into bins and merging adjacent regions that fall into the same bin even if they have
different exact coverage values. This can dramatically reduce the size of the output compared to the per-base.

It also allows outputting regions of low, high, and "callable" coverage as in [GATK's callable loci tool](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_coverage_CallableLoci.php).

An example of quantize arguments:
```
--quantize 0:1:4:100:200: # ... arbitary number of quantize bins.
```

indicates bins of: 0, 1-3, 4-99, 100-200, 200-infinity

The default for `mosdepth` is to output the integer index of the bin associated with each interval. So in the example here,
a depth of 0 would be assigned to bin 0, a depth of 3 to bin 1, a depth of 5 to bin 2 and so-on.

To change what is reported as the bin number, a user can set environment variables e.g.:

```
export MOSDEPTH_Q0=NO_COVERAGE
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=CALLABLE
export MOSDEPTH_Q3=HIGH_COVERAGE
```

In this case, the bin number is replaced by the text in the appropriate environment variable.

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
length, so for the 249MB chromosome 1, it will require 1GB of memory.

`mosdepth` is written in [nim](https://nim-lang.org/) and it uses our [htslib](https://github.com/samtools/htslib)
via our nim wrapper [hts-nim](https://github.com/brentp/hts-nim/)

## speed and memory comparison

`mosdepth`, `samtools`, `bedtools`, and `sambamba` were run on a 30X genome.
relative times are relative to mosdepth per-base mode with a single thread.

`mosdepth` can report the mean depth in 500-base windows genome-wide info
under 9 minutes of user time with 3 threads.

| format |    tool    | threads  | mode   | relative time | run-time | memory |
| ------ | ---------- | -------- | ------ | ------------- | -------  | -------|
|  BAM   |  mosdepth  |    1     | base   |     1         |  25:23   |  1196  |
|  BAM   |  mosdepth  |    3     | base   |    0.57       |  14:27   |  1197  |
|  CRAM  |  mosdepth  |    1     | base   |    1.17       |  29:47   |  1205  |
|  CRAM  |  mosdepth  |    3     | base   |    0.56       |  14:08   |  1225  |
|  BAM   |  mosdepth  |    3     | window |    0.34       |  8:44    |  1277  |
|  BAM   |  mosdepth  |    1     | window |    0.80       |  20:26   |  1212  |
|  CRAM  |  mosdepth  |    3     | window |    0.35       |  8:47    |  1233  |
|  CRAM  |  mosdepth  |    1     | window |    0.88       |  22:23   |  1209  |
|  BAM   |  sambamba  |    1     | base   |    5.71       | 2:24:53  |  166   |
|  BAM   |  samtools  |    1     | base   |    1.98       | 50:12    |  27    |
|  CRAM  |  samtools  |    1     | base   |    1.79       | 45:21    |  451   |
|  BAM   |  bedtools  |    1     | base   |    5.31       | 2:14:44  |  1908  |


Note that the threads to `mosdepth` (and samtools) are decompression threads. After
about 4 threads, there is no benefit for additional threads:

![mosdepth-scaling](https://user-images.githubusercontent.com/1739/31246294-256d1b7c-a9ca-11e7-8e28-6c4d07cba3f5.png)


### Accuracy

We compared `samtools depth` with default arguments to `mosdepth` without overlap detection and discovered **no
differences across the entire chromosome**.
