# mosdepth

fast depth calculation and normalization

`mosdepth` can output per-base depth about twice as fast `samtools depth`.

it can also output per-window depth given a window size or a BED file of regions.

## usage

```
nimdepth

  Usage: nimdepth [options] <BAM-or-CRAM>

  -t --threads <threads>     number of CRAM|BAM decompression threads to use [default: 0]
  -c --chrom <chrom>         chromosome to restrict depth calculation.
  -Q --mapq <mapq>           mapping quality threshold [default: 0]
  -b --by <bed|window>       BED file of regions or an (integer) window-size for which to calculate depth.
  -f --fasta <fasta>         fasta file for use with CRAM files.
  -h --help                  show help
```

### exome example

To calculate the coverage in each exome capture region:
```
mosdepth -b capture.bed sample.exome.bam > sample.exome.coverage.bed
```
For a 5.5GB exome file and 1,195,764 capture regions, this completes
in XXX with a single CPU.

### WGS example

For per-base whole-genome coverage:

```
mosdepth sample.wgs.bam > sample.txt
```

For 500-base windows:

```
mosdepth --by 500 sample.wgs.bam > sample.500.bed
```

## how it works

For each chromosome, `mosdepth` creates an array the length of the chromosome. For every
start it encounters, integer increments value in that position of the array. For every
stop, it decrements that position. From this, the depth at a particular position is the
cumulative sum of all array positions preceding it (a similar algorithm is used in BEDTools
where starts and stops are tracked separately). `mosdepth` avoids double-counting
overlapping mate-pairs and it tracks every aligned part of every read using the CIGAR
oparations.

This array accounting is very fast. There are no extra allocations or objects to track and
it is also conceptually simple. For these regions, it is faster than `samtools depth` which
works by using the [pileup](http://samtools.sourceforge.net/pileup.shtml) machinery that
tracks each read, each base. 

This method has some limitations. Because a large array is allocated and we (in general)
must take the cumulative sum of all preceding positions to know the depth at any position.
As such, it is slower for small, 1-time regional queries. It is, however fast for window-based
or BED-based regions, because it first calculates the full chromosome coverage and then
reports the coverage for each region in that chromosome.
