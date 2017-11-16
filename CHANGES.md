0.2.1
=====
+ allow unsorted bed as input to --by

0.2.0
=====
+ **2X speed improvement for CRAM** by not decoding unused base-qualities.
+ add new `--thresholds` argument. See README for usage. thresholds and quantization are highly recommended over
+ when using quantize, labels now indicate the range of depths encompassed by that region.
+ for quantize and per-base, mosdepth will output all chromosomes in the bam header, even if they have no alignments.
  this makes it so that mosdepth output for files aligned to the same reference have the same number and order of chromosomes
  and same total base coverage
  per-base output when possible as they are more compact and faster to output.
+ fix bug with dist output for chroms larger than 2.1 billion bases that also affected total output in some cases.
