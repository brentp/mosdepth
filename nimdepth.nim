import hts as hts
import tables
import strutils as S
import algorithm as alg
import sequtils as sequtils
import strutils as su
import os
import docopt
import tables

proc setvbuf(stream: File, buf: cstring, buftype: int, size: int32): int {.importc: "setvbuf", header:"<stdio.h>".}
type
  pair = tuple[pos: int, value: int32]
  depth_t = tuple[pos: int, value: int]
  region_t = ref object
    chrom: string
    start: int
    stop: int

proc length(r: region_t): int =
  return r.stop - r.start

proc `$`(r: region_t): string =
  if r.stop != 0:
    return format("$1:$2-$3", r.chrom, r.start + 1, r.stop)
  else:
    return format("$1:$2", r.chrom, r.start + 1)

iterator gen_depths(arr: seq[int32], window: int, stop: var int): depth_t =
  # given `arr` with values in each index indicating the number of reads
  # starting or ending at that location, generate depths.
  var
    last_start = -1
    last_depth = -1
    depth = 0
    i = -1
  if stop <= 0:
    stop = len(arr)
  
  for change in arr:
    inc(i)
    depth += int(change)
    if i == stop:
      break
    if depth == last_depth and (window == 0 or (i mod window) != 0): continue
    if last_depth != -1:
      yield (i, last_depth)

    last_start = i
    last_depth = depth
  if last_depth != -1:
    yield (i, last_depth)

proc scaled_sum(depths: seq[int], weights: var seq[float64], n: int, window: float64): float64 =
  for i, x in weights[0 .. <n]:
    weights[i] = x / window
  result = float64(0)
  var i = 0
  for d in depths[0 .. < n]:
    result += float64(d) * weights[i]
    inc(i)

proc by_window(arr: seq[int32], chrom: string, window: int, region: region_t) =
  # discretize the depth by window.
  # in some cases, window will split a region and in other cases, a region
  # will define the window. e.g. region = chr1:1000-2500 and window = 1500
  var weights = new_seq[float64](window)
  var depths = new_seq[int](window)
  var stop = len(arr)
  var start = 0
  if region != nil:
    stop = region.stop
    start = region.start
  var last_pos = 0
  var n = 1
  var tchrom = chrom & "\t"
  for p in gen_depths(arr, window, stop):
    depths[n] = p.value
    weights[n] = float64(p.pos - last_pos)
    if (p.pos - start) mod window == 0:
      var dp = su.format_float(scaled_sum(depths, weights, n, float64(window)), ffDecimal, precision=3)
      su.trim_zeros(dp)
      stdout.write_line tchrom & intToStr(p.pos - window) & "\t" & intToStr(p.pos) & "\t" & dp
      n = 0
      if stop != 0 and p.pos == stop:
        return
    inc(n)
    if region != nil and region.stop != 0 and p.pos >= region.stop:
      break
    last_pos = p.pos
  if n > 0:
    var left = (stop mod window)
    var dp = su.format_float(scaled_sum(depths, weights, n, float64(left)), ffDecimal, precision=3)
    su.trim_zeros(dp)
    stdout.write_line tchrom & intToStr(stop - left) & "\t" & intToStr(stop) & "\t" & dp

proc dump(arr: seq[int32], chrom: string, region: region_t) =
  var tchrom = chrom & "\t"
  var stop = 0
  if region != nil and region.stop != 0:
    stop = region.stop
  for p in gen_depths(arr, 0, stop):
    stdout.write_line tchrom & intToStr(p.pos) & "\t" & intToStr(p.value)

proc write_depth(arr: seq[int32], chrom: string, window: int, region: region_t) =
  if window == 0:
    dump(arr, chrom, region)
  else:
    by_window(arr, chrom, window, region)

proc pair_sort(a, b: pair): int =
   return a.pos - b.pos

iterator gen_start_ends(c: Cigar, ipos: int): pair =
  # generate start, end pairs given a cigar string and a position offset.
  if c.len == 1 and c[0].op == CigarOp(match):
    yield (ipos, int32(1))
    yield (ipos + c[0].len, int32(-1))
  else:
    var pos = ipos
    var last_stop = 0
    var con: Consume
    for op in c:
      con = op.consumes
      if not con.reference:
        continue
      var olen = op.len
      if con.query:
        if pos != last_stop:
          yield (pos, int32(1))
          if last_stop != 0:
            yield (last_stop, int32(-1))
        last_stop = pos + olen
      pos += olen
    if last_stop != 0:
      yield (last_stop, int32(-1))

proc inc_coverage(c: Cigar, ipos: int = 0, arr: var seq[int32]) {. inline .} =
  for p in gen_start_ends(c, ipos):
    arr[p.pos] += p.value

iterator regions(bam: hts.Bam, region: region_t): Record =
  if region == nil:
    for r in bam:
      yield r
  elif region != nil:
    for r in bam.querys($region):
      yield r

proc bed_line_to_region(line: string): region_t =
   var cse = sequtils.to_seq(line.strip().split("\t"))
   var s = S.parse_int(cse[1])
   var e = S.parse_int(cse[2])
   var reg = region_t(chrom: cse[0], start:s, stop: e)
   return reg

proc region_line_to_region(region: string): region_t =
  if region == nil:
    return nil
  var i = 0
  var r = region_t()
  for w in region.split({':', '-'}):
    if i == 1:
      r.start = S.parse_int(w) - 1
    elif i == 2:
      r.stop = S.parse_int(w)
    else:
      r.chrom = w
    inc(i)
  return r

proc main(bam: hts.Bam, arr: var seq[int32], region: region_t, mapq:int= -1, window: int=0) =
  var seqs = bam.hdr.targets

  var tgt: hts.Target
  var mate: Record
  var seen = newTable[string, Record]()

  for rec in bam.regions(region):
    if rec.flag.has_flag(BAM_FUNMAP or BAM_FQCFAIL or BAM_FDUP or BAM_FSECONDARY): continue
    if int(rec.qual) < mapq: continue
    if tgt == nil or tgt.tid != rec.b.core.tid:
        if tgt != nil:
          write_depth(arr, tgt.name, window, region)
        tgt = seqs[rec.b.core.tid]
        arr = new_seq[int32](tgt.length)
        seen.clear()
    # rec:   --------------
    # mate:             ------------
    # handle overlapping mate pairs.
    if rec.flag.proper_pair:
      if rec.stop > rec.matepos and rec.start < rec.matepos:
        seen[rec.qname] = rec.copy()
      else:
        if seen.take(rec.qname, mate):
          # we have an overlapping pair, and we know that mate is lower. e.g
          # mate:   --------------
          # rec:             ------------
          # decrement:       -----
          if rec.b.core.n_cigar == 1 and mate.b.core.n_cigar == 1:
            dec(arr[rec.start])
            inc(arr[mate.stop])
          else:
            # track the overlaps of pair.
            # anywhere there is overlap, the cumulative sum of pair.depth will be 2. we dec the start and inc the end of the overlap.
            # this removes the double counting.
            # e.g.:
            # (pos: 4623171, value: 1)(pos: 4623223, value: 1)(pos: 4623240, value: -1)(pos: 4623241, value: 1)(pos: 4623264, value: -1)(pos: 4623320, value: -1)
            # would dec the following intervals:
            # 4623223 4623240
            # 4623241 4623264
            # chr1 4623171 69M1D23M9S (pos: 4623171, value: 1)(pos: 4623241, value: 1)(pos: 4623240, value: -1)(pos: 4623264, value: -1)
            # chr1 4623223 4S97M (pos: 4623223, value: 1)(pos: 4623320, value: -1)
            assert rec.start < mate.stop
            # each element will have a .value of 1 for start and -1 for end.

            var ses = sequtils.to_seq(gen_start_ends(rec.cigar, rec.start))
            for p in gen_start_ends(mate.cigar, mate.start):
                ses.add(p)
            alg.sort(ses, pair_sort)
            var pair_depth = 0
            var last_pos = 0
            #if len(ses) > 4: stderr.write_line ses
            for p in ses:
              assert pair_depth <= 2
              # if we are at pair_depth 2, there is overlap and when the incoming
              # value is -1, then it is dropping back down to 1.
              if p.value == -1 and pair_depth == 2:
                #if len(ses) > 4: stderr.write_line last_pos, " ", p.pos
                dec(arr[last_pos])
                inc(arr[p.pos])
              pair_depth += p.value
              last_pos = p.pos
            assert pair_depth == 0

    inc_coverage(rec.cigar, rec.start, arr)
  if tgt != nil:
    write_depth(arr, tgt.name, window, region)
  elif region != nil:
    write_depth(arr, region.chrom, window, region)

proc bed_main(bam: hts.Bam, bed: string, thread:int=0, mapq:int= -1, window:int=0) =
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  var arr: seq[int32]
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    var rl = bed_line_to_region($kstr.s)
    main(bam, arr, rl, mapq, window=rl.length)

  hts.free(kstr.s)

when(isMainModule):

  let doc = """
  nimdepth

  Usage: nimdepth [options] <BAM>
  
  -t --threads <threads>  number of threads to use [default: 0]
  -Q --mapq <mapq>        mapping quality threshold [default: 0]
  -r --region <region>    region to restrict depth calc.
  -L --bed <bed>          BED file of regions for which to calculate depth.
                          Can not be used with `--region`.
  -w --window <window>    window-size.
  -h --help               show help
  """

  let args = docopt(doc, version = "nimdepth 0.1.1")
  let mapq = S.parse_int($args["--mapq"])
  var window = 0
  var bed_based = false 
  if $args["--window"] != "nil":
    window = S.parse_int($args["--window"])
  var region: string
  if $args["--region"] != "nil":
    region = $args["--region"]
    assert $args["--bed"] == "nil"
  elif $args["--bed"] != "nil":
    bed_based = true
  GC_disableMarkAndSweep()
  discard setvbuf(stdout, nil, 0, 16384)
  var arr : seq[int32]

  if bed_based:
    var bam = hts.open_hts($args["<BAM>"], threads=1, index=true)
    bed_main(bam, $args["--bed"], S.parse_int($args["--threads"]), mapq, window=window)
  else:
    var threads = S.parse_int($args["--threads"])
    var bam = hts.open_hts($args["<BAM>"], threads=threads, index=region != nil)
    main(bam, arr, region_line_to_region(region), mapq, window=window)
