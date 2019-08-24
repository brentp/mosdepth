import hts
import tables
import strutils as S
import algorithm as alg
import sequtils as sequtils
import strutils as su
import os
import docopt
import times
import math
import ./depthstat

var precision: int
var output_summary_header = true

try:
  var tmp = getEnv("MOSDEPTH_PRECISION")
  precision = parse_int(tmp)
except:
  precision = 2

type
  pair = tuple[pos: int, value: int32]
  depth_t = tuple[start: int, stop: int, value: int]
  depth_s = tuple[start: int, stop: int, value: string]
  region_t = ref object
    chrom: string
    start: uint32
    stop: uint32
    name: string

  coverage_t = seq[int32]

proc `$`*(r: region_t): string =
  if r == nil:
    return ""
  if r.stop != 0:
    return format("$1:$2-$3", r.chrom, r.start + 1, r.stop)
  else:
    return format("$1:$2", r.chrom, r.start + 1)

proc to_coverage(c: var coverage_t) =
  # to_coverage converts from an array of start/end inc/decs to actual coverage.
  var cum = int32(0)
  for i, d in c.mpairs:
    cum += d
    d = cum

iterator gen_depths(arr: coverage_t, offset: int=0, istop: int=0): depth_t =
  # given `arr` with values in each index indicating the number of reads
  # starting or ending at that location, generate depths.
  # offset is only used for a region like chr6:200-30000, in which case, offset will be 200
  var
    last_depth = -1
    i = 0
    last_i = 0
    stop: int
  if istop <= 0:
    stop = len(arr)-1
  else:
    stop = istop
  # even with an offset, have to start from the beginning of the array
  # to get the proper depth.
  for depth in arr:
    if i == stop:
      break
    if i < offset or depth == last_depth:
      inc(i)
      continue

    if last_depth != -1:
      yield (last_i, i, last_depth)

    last_depth = depth
    last_i = i
    if i + 1 == stop: break
    inc(i)

  if last_i < stop:
    yield (last_i, len(arr)-1, last_depth)

  # this is horrible, but it works. we don't know
  # if we've already printed the record on not.
  elif  last_i != i:
      yield (last_i - 1, i, last_depth)
  else:
      yield (last_i, i, last_depth)

proc linear_search*(q:int, vals:seq[int], idx: ptr int) {.inline.} =
  if q < vals[0] or q > vals[vals.high]:
    idx[] = -1
    return
  for i, val in vals:
    if val > q:
      idx[] = i - 1
      return
    if val == q:
      idx[] = i
      return
  idx[] = vals.high

proc make_lookup*(quants: seq[int]): seq[string] =
  var L = new_seq[string](len(quants)-1)
  for i in 0..L.high:
    var t = getEnv("MOSDEPTH_Q" & intToStr(i))
    if t == "":
      if quants[i+1] == high(int):
        L[i] = intToStr(quants[i]) & ":inf"
      else:
        L[i] = intToStr(quants[i]) & ":" & intToStr(quants[i+1])
    else:
      L[i] = t
  return L

iterator gen_quantized(quants: seq[int], arr: coverage_t): depth_s {.inline.} =
  # like gen_depths but merges adjacent entries in same quantize bins.
  if len(arr) > 0:
    var lookup = make_lookup(quants)
    var last_quantized, quantized: int
    linear_search(arr[0], quants, last_quantized.addr)
    var last_pos = 0
    # slicing into the array does a copy.
    for pos in 0..<(arr.high-1):
      let d = arr[pos]
      linear_search(d, quants, quantized.addr)
      if quantized == last_quantized: continue
      if last_quantized != -1 and last_quantized < len(lookup):
        yield (last_pos, pos, lookup[last_quantized])
      last_quantized = quantized
      last_pos = pos
    if last_quantized != -1 and last_pos < arr.high and last_quantized < len(lookup):
      yield (last_pos, len(arr)-1, lookup[last_quantized])

proc pair_sort(a, b: pair): int =
   return a.pos - b.pos

iterator gen_start_ends(c: Cigar, ipos: int): pair {.inline.} =
  # generate start, end pairs given a cigar string and a position offset.
  if c.len == 1 and c[0].op == CigarOp.match:
    yield (ipos, int32(1))
    yield (ipos + c[0].len, int32(-1))
  else:
    var pos = ipos
    var last_stop = -1
    var con: Consume
    for op in c:
      con = op.consumes
      if not con.reference:
        continue
      var olen = op.len
      if con.query:
        if pos != last_stop:
          yield (pos, int32(1))
          if last_stop != -1:
            yield (last_stop, int32(-1))
        last_stop = pos + olen
      pos += olen
    if last_stop != -1:
      yield (last_stop, int32(-1))

proc inc_coverage(c: Cigar, ipos: int = 0, arr: var seq[int32]) {.inline.} =
  for p in gen_start_ends(c, ipos):
      arr[p.pos] += p.value

iterator regions(bam: hts.Bam, region: region_t, tid: int, targets: seq[hts.Target]): Record {.inline.} =
  if region == nil:
    for r in bam:
      yield r
  elif region != nil:
    var stop = region.stop
    if tid != -1:
      if stop == 0:
        stop = targets[tid].length
      for r in bam.query(tid, int(region.start), int(stop)):
        yield r
    else:
      stderr.write_line("[mosdepth]", region.chrom, " not found")

proc bed_line_to_region(line: string): region_t =
   var
     cse = line.strip().split('\t', 5)

   if len(cse) < 3:
     stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
     return nil
   var
     s = S.parse_int(cse[1])
     e = S.parse_int(cse[2])
     reg = region_t(chrom: cse[0], start: uint32(s), stop: uint32(e))
   doAssert s <= e, "[slivar] ERROR: start > end in bed line:" & line
   if len(cse) > 3:
     reg.name = cse[3]
   return reg

proc region_line_to_region(region: string): region_t =
  if region == "" or region == "nil":
    return nil
  var i = 0
  var r = region_t()
  for w in region.split({':', '-'}):
    if i == 1:
      r.start = uint32(S.parse_int(w)) - 1
    elif i == 2:
      r.stop = uint32(S.parse_int(w))
    else:
      r.chrom = w
    inc(i)
  return r

proc get_tid(tgts: seq[hts.Target], chrom: string): int =
  for t in tgts:
    if t.name == chrom:
      return t.tid

proc init(arr: var coverage_t, tlen:int) =
  ## try to re-use the array.
  if len(arr) != int(tlen):
    # must create a new array in some cases.
    if arr.len == 0:
      arr = new_seq[int32](tlen)
      return
    else:
      # otherwise can re-use and zero
      arr.set_len(int(tlen))
  zeroMem(arr[0].addr, len(arr) * sizeof(arr[0]))

proc coverage(bam: hts.Bam, arr: var coverage_t, region: var region_t, mapq:int= -1, eflag: uint16=1796, iflag:uint16=0, read_groups:seq[cstring]=(@[]), fast_mode:bool=false): int =
  # depth updates arr in-place and yields the tid for each chrom.
  # returns -1 if the chrom is not found in the bam header
  # returns -2 if the chrom was found in the header, but there was no data for it
  # otherwise returns the tid.
  var
    targets = bam.hdr.targets
    tgt: hts.Target
    mate: Record
    seen = newTable[string, Record]()
    has_read_groups = read_groups.len > 0

  var tid = if region != nil: get_tid(targets, region.chrom) else: -1
  if tid == -1:
    return -1

  tgt = targets[tid]

  var found = false
  for rec in bam.regions(region, tid, targets):
    if not found:
      arr.init(int(tgt.length+1))
      found = true
    if int(rec.mapping_quality) < mapq: continue
    if (rec.flag and eflag) != 0:
      continue
    if iflag != 0 and ((rec.flag and iflag) == 0):
      continue
    if has_read_groups:
      var t = tag[cstring](rec, "RG")
      if t.isNone or not read_groups.contains(t.get):
        continue
    if tgt.tid != rec.b.core.tid:
        raise newException(OSError, "expected only a single chromosome per query")

    # rec:   --------------
    # mate:             ------------
    # handle overlapping mate pairs.
    if (not fast_mode) and rec.flag.proper_pair and (not rec.flag.supplementary):
      if rec.b.core.tid == rec.b.core.mtid and rec.stop > rec.matepos and 
        # First case is partial overlap, second case is complete overlap
        # For complete overlap we must check if the mate was already seen or not yet
        ((rec.start < rec.matepos) or (rec.start == rec.mate_pos and not seen.hasKey(rec.qname))):
        var rc = rec.copy()
        seen[rc.qname] = rc
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
            assert (rec.start <= mate.stop), rec.tostring() & "\n" & mate.tostring()
            # each element will have a .value of 1 for start and -1 for end.

            var ses = sequtils.to_seq(gen_start_ends(rec.cigar, rec.start))
            for p in gen_start_ends(mate.cigar, mate.start):
                ses.add(p)
            alg.sort(ses, pair_sort)
            var pair_depth = 0
            var last_pos = 0
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
            if pair_depth != 0: echo $rec.qname & ":" & $rec & " " & $mate.qname & ":" & $mate & " " & $pair_depth
    if fast_mode:
      arr[rec.start] += 1
      arr[rec.stop] -= 1
    else:
      inc_coverage(rec.cigar, rec.start, arr)

  if not found:
    return -2
  return tgt.tid

proc bed_to_table(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var kstr = kstring_t(l:0, m: 0, s: nil)
  var hf = hts_open(cstring(bed), "r")
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if kstr.s[0] == 't' and ($kstr.s).startswith("track "):
      continue
    if kstr.s[0] == '#':
      continue
    var v = bed_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
      sort(ivs, proc (a, b: region_t): int = int(a.start) - int(b.start))

  hts.free(kstr.s)
  return bed_regions

iterator window_gen(window: uint32, t: hts.Target): region_t =
  var start:uint32 = 0
  while start + window < t.length:
    yield region_t(chrom: t.name, start: start, stop: start + window)
    start += window
  if start != t.length:
    yield region_t(chrom: t.name, start: start, stop: t.length)

iterator region_gen(window: uint32, target: hts.Target, bed_regions: TableRef[string, seq[region_t]]): region_t =
    if bed_regions == nil:
      for r in window_gen(window, target): yield r
    else:
      if bed_regions.contains(target.name):
        for r in bed_regions[target.name]: yield r
        bed_regions.del(target.name)

proc imean(vals: coverage_t, start:uint32, stop:uint32, ms:var CountStat[uint32]): float64 =
  if start > uint32(len(vals)):
    return 0

  if ms.len != 0:
    ms.clear()
    for i in start..<min(stop, uint32(len(vals))):
      ms.add(vals[i])
    return ms.median.float64

  else:
    var L = float64(stop - start)
    for i in start..<min(stop, uint32(len(vals))):
      result += float64(vals[int(i)]) / L

const MAX_COVERAGE = int32(400000)

proc inc(d: var seq[int64], coverage: var coverage_t, start:uint32, stop:uint32) =
  var v:int32
  var L = int32(d.high)
  if int(start) >= len(coverage):
    stderr.write_line("[mosdepth] warning requested interval outside of chromosome range:", start, "..", stop)
    return
  var istop = min(stop, uint32(coverage.len))

  for i in start..<istop:
    v = coverage[i]
    if v > MAX_COVERAGE:
      v = MAX_COVERAGE - 10
    if v >= L:
      d.set_len(v + 10)
      for j in (L+1).int..d.high:
        d[j] = 0
      L = int32(d.high)
    if v < 0: continue
    d[v] += 1

proc write_distribution(chrom: string, d: var seq[int64], fh:File) =
  var sum: int64
  for v in d:
    sum += int64(v)

  if sum < 1: return
  var cum: float64 = 0
  # reverse and then cumsum so that e.g. a value of 1 is the proportion of
  # bases with a coverage of at least 1.
  reverse(d)
  # skip until we see a non-zero value for high-coverage end of array.
  for i, v in pairs(d):
    var irev = len(d) - i - 1
    if irev > 300 and v == 0: continue
    cum += float64(v) / float64(sum)
    if cum < 8e-5: continue
    fh.write_line(chrom, "\t", $irev & "\t" & su.format_float(cum, ffDecimal, precision=precision))
  # reverse it back because we use to update the full genome
  reverse(d)

proc write_summary(region: string, stat: depth_stat, fh:File) =
  var mean_depth: float64
  if stat.cum_length > 0:
    mean_depth = float64(stat.cum_depth) / float64(stat.cum_length)
  else:
    mean_depth = 0.float64
  let stat_min = if stat.min_depth == uint32.high: 0.uint32 else: stat.min_depth
  if output_summary_header:
    fh.write_line ["chrom",
                   "length",
                   "bases",
                   "mean",
                   "min",
                   "max"].join("\t")
    output_summary_header = false
  fh.write_line [region,
                 $stat.cum_length,
                 $stat.cum_depth,
                 $mean_depth.format_float(ffDecimal, precision=precision),
                 $stat_min,
                 $stat.max_depth].join("\t")

proc get_targets(targets: seq[hts.Target], r: region_t): seq[hts.Target] =
  if r == nil:
    return targets
  result = new_seq[hts.Target](1)
  for t in targets:
    if t.name == r.chrom:
      result[0] = t

# copy the chromosome array into the genome-wide
proc sum_into(afrom: seq[int64], ato: var seq[int64]) =
  if len(afrom) > len(ato):
    var b = len(ato)
    ato.set_len(afrom.len)
    for i in b..ato.high:
      ato[i] = 0

  for i in 0..afrom.high:
    ato[i] += afrom[i]

proc get_quantize_args*(qa: string) : seq[int] =
  if qa == "nil":
    return
  var a = qa
  if a.count(':') == 0:
    a = ':' & a & ':'
  # if it starts with : we go prepend 0
  if a[0] == ':':
    a = "0" & a
  # if it ends with ":" we make that bin include all numbers above it.
  if a[a.high] == ':':
    a = a & intToStr(high(int))
  try:
    var qs = map(a.split(':'), proc (s:string): int = return parse_int(s))
    sort(qs, system.cmp[int])
    return qs
  except:
    stderr.write_line("[mosdepth] invalid quantize string: '" & a & "'")
    quit(2)


proc write_thresholds(fh:BGZI, tid:int, arr:var coverage_t, thresholds:seq[int], region: region_t) =
  # write the number of bases in each region that are >= each threshold.
  if thresholds.len == 0: return
  var
    line = new_string_of_cap(32)
    start = int(region.start)
    stop = int(region.stop)
  line.add(region.chrom & "\t")
  line.add(intToStr(start) & "\t")
  line.add(intToStr(stop))
  if region.name != "":
    line.add("\t" & region.name)
  else:
    line.add("\tunknown")

  if tid == -2:
    for i in thresholds:
      line.add("\t0")
    discard fh.write_interval(line, region.chrom, start, stop)
    return

  var counts = new_seq[int](len(thresholds))
  shallow(arr)

  # iterate over the region and count bases >= request cutoffs.
  for v in arr[start..<stop]:
    for i, t in thresholds:
      # if we know they are sorted we can break
      if v < t: break
      counts[i] += 1

  for count in counts:
    line.add("\t" & intToStr(count))
  discard fh.write_interval(line, region.chrom, start, stop)

proc write_header(fh:BGZI, thresholds: seq[int]) =
  discard fh.bgz.write("#chrom	start	end	region")
  for threshold in thresholds:
    discard fh.bgz.write("\t" & intToStr(threshold) & "X")
  discard fh.bgz.write("\n")

proc get_min_levels(targets: seq[Target]): int =
  # determine how many levels are needed to store the data given
  # the largest chromosome
  var max_len = targets[0].length.uint64
  for t in targets:
    if t.length > max_len:
      max_len = t.length.uint64

  result = 0
  var s = (1 shl 14).uint64
  while max_len > s:
    result += 1
    s = s shl 3


proc main(bam: hts.Bam, chrom: region_t, mapq: int, eflag: uint16, iflag: uint16, region: string, thresholds: seq[int],
          fast_mode:bool, args: Table[string, docopt.Value], use_median:bool=false) =
  # windows are either from regions, or fixed-length windows.
  # we assume the input is sorted by chrom.
  var
    targets = bam.hdr.targets
    sub_targets = get_targets(targets, chrom)
    read_groups: seq[cstring]
    rchrom : region_t
    arr: coverage_t
    prefix: string = $(args["<prefix>"])
    skip_per_base = args["--no-per-base"]
    window: uint32 = 0
    bed_regions: TableRef[string, seq[region_t]] # = Table[string, seq[region_t]]
    fbase: BGZI
    #fbase: BGZ
    fquantize: BGZI
    fthresholds: BGZI
    fregion: BGZI
    fh_global_dist:File
    fh_region_dist:File
    fh_summary: File
    quantize = get_quantize_args($args["--quantize"])

    # summary stat output
    chrom_region_stat: depth_stat
    chrom_stat: depth_stat
    global_region_stat: depth_stat
    global_stat: depth_stat

  var region_distribution = new_seq[int64](1000)
  var global_distribution = new_seq[int64](1000)

  if $args["--read-groups"] != "nil":
    for r in ($args["--read-groups"]).split(','):
      read_groups.add(r.cstring)
  var levels = get_min_levels(targets)

  var chrom_region_distribution, chrom_global_distribution: seq[int64]

  if not skip_per_base:
    # can't use set-threads when indexing on the fly so this must
    # not call set_threads().
    fbase = wopen_bgzi(prefix & ".per-base.bed.gz", 1, 2, 3, true, compression_level=1, levels=levels)
    #open(fbase, prefix & ".per-base.bed.gz", "w1")
  if quantize.len != 0:
    fquantize = wopen_bgzi(prefix & ".quantized.bed.gz", 1, 2, 3, true, compression_level=1, levels=levels)

  if thresholds.len != 0:
    fthresholds = wopen_bgzi(prefix & ".thresholds.bed.gz", 1, 2, 3, true, compression_level=1, levels=levels)
    fthresholds.write_header(thresholds)

  if not open(fh_global_dist, prefix & ".mosdepth.global.dist.txt", fmWrite):
    stderr.write_line("[mosdepth] could not open file:", prefix & ".mosdepth.global.dist.txt")

  if not open(fh_summary, prefix & ".mosdepth.summary.txt", fmWrite):
    stderr.write_line("[mosdepth] could not open file:", prefix & ".mosdepth.summary.txt")

  if region != "" and not open(fh_region_dist, prefix & ".mosdepth.region.dist.txt", fmWrite):
    stderr.write_line("[mosdepth] could not open file:", prefix & ".mosdepth.dist.txt")

  if region != "":
    fregion = wopen_bgzi(prefix & ".regions.bed.gz", 1, 2, 3, true, levels=levels)
    if region.isdigit():
      window = uint32(S.parse_int(region))
    else:
      bed_regions = bed_to_table(region)
  shallow(arr)

  var cs = initCountStat[uint32](size=if use_median: 65536 else: 0)

  for target in sub_targets:
    chrom_global_distribution = new_seq[int64](1000)
    if region != "":
      chrom_region_distribution = new_seq[int64](1000)
    # if we can skip per base and there's no regions from this chrom we can avoid coverage calc.
    if skip_per_base and thresholds.len == 0 and quantize.len == 0 and bed_regions != nil and not bed_regions.contains(target.name):
      continue
    rchrom = region_t(chrom: target.name)
    var tid = coverage(bam, arr, rchrom, mapq, eflag, iflag, read_groups=read_groups, fast_mode=fast_mode)
    if tid == -1: continue # -1 means that chrom is not even in the bam
    if tid != -2: # -2 means there were no reads in the bam
      arr.to_coverage()

    var starget = target.name & "\t"
    if region != "":
      var line = new_string_of_cap(16384)
      var me = 0'f64
      for r in region_gen(window, target, bed_regions):
        if tid != -2:
          me = imean(arr, r.start, r.stop, cs)
          chrom_region_stat = chrom_region_stat + newDepthStat(arr[r.start..<r.stop])
        var m = su.format_float(me, ffDecimal, precision=precision)

        if r.name == "":
          line.add(starget & intToStr(int(r.start)) & "\t" & intToStr(int(r.stop)) & "\t" & m)
        else:
          line.add(starget & intToStr(int(r.start)) & "\t" & intToStr(int(r.stop)) & "\t" & r.name & "\t" & m)
        discard fregion.write_interval(line, target.name, int(r.start), int(r.stop))
        line = line[0..<0]
        if tid != -2:
          chrom_region_distribution.inc(arr, r.start, r.stop)
        write_thresholds(fthresholds, tid, arr, thresholds, r)
    if tid != -2:
      chrom_global_distribution.inc(arr, uint32(0), uint32(len(arr) - 1))
      chrom_stat = newDepthStat(arr[0..<len(arr)-1])
      global_stat = global_stat + chrom_stat
      write_summary(target.name, chrom_stat, fh_summary)
      if region != "":
        write_summary(target.name & "_region", chrom_region_stat, fh_summary)
      global_region_stat = global_region_stat + chrom_region_stat
      chrom_region_stat.clear()

    # write the distribution for each chrom
    write_distribution(target.name, chrom_global_distribution, fh_global_dist)
    if region != "":
      write_distribution(target.name, chrom_region_distribution, fh_region_dist)

    # then copy it to the genome distribution
    if tid >= 0:
      if region != "":
        sum_into(chrom_region_distribution, region_distribution)
      sum_into(chrom_global_distribution, global_distribution)

    if not skip_per_base:
      if tid == -2:
        discard fbase.write_interval(starget & "0\t" & intToStr(int(target.length)) & "\t0", target.name, 0, int(target.length))
      else:
        for p in gen_depths(arr):
          discard fbase.write_interval(starget & intToStr(p.start) & "\t" & intToStr(p.stop) & "\t" & intToStr(p.value), target.name, p.start, p.stop)
    if quantize.len != 0:
      if tid == -2 and quantize[0] == 0:
        var lookup = make_lookup(quantize)
        discard fquantize.write_interval(starget & "0\t" & intToStr(int(target.length)) & "\t" & lookup[0], target.name, 0, int(target.length))
      else:
        if tid == -2: continue
        for p in gen_quantized(quantize, arr):
            discard fquantize.write_interval(starget & intToStr(p.start) & "\t" & intToStr(p.stop) & "\t" & p.value, target.name, p.start, p.stop)

  write_summary("total", global_stat, fh_summary)
  if region != "":
    write_summary("total_region", global_region_stat, fh_summary)

  write_distribution("total", global_distribution, fh_global_dist)
  if region != "":
    write_distribution("total", region_distribution, fh_region_dist)
    fh_region_dist.close()

  if bed_regions != nil and chrom == nil:
    for chrom, regions in bed_regions:
      stderr.write_line("[mosdepth] warning chromosome:", chrom, " from bed with " , len(regions), " regions not found")

  if fregion != nil and close(fregion) != 0:
      stderr.write_line("[mosdepth] error writing region file\n")
      quit(1)

  if fquantize != nil and close(fquantize) != 0:
      stderr.write_line("[mosdepth] error writing quantize file\n")
      quit(1)

  if fthresholds != nil and close(fthresholds) != 0:
      stderr.write_line("[mosdepth] error writing thresholds file\n")
      quit(1)

  if fbase != nil and close(fbase) != 0:
      stderr.write_line("[mosdepth] error writing per-base file\n")
      quit(1)
  close(fh_global_dist)

proc check_chrom(r: region_t, targets: seq[Target]) =
  if r == nil: return
  for t in targets:
    if t.name == r.chrom:
      return
  stderr.write_line "[mosdepth] chromosome ", r.chrom, " not found"
  quit(1)

proc threshold_args*(ts: string): seq[int] =
  if ts == "nil":
    return
  result = map(ts.split(','), proc (s:string): int = return parse_int(s))
  sort(result)


proc check_cram_has_ref(cram_path: string, fasta:string) =
  if fasta != "" and exists_file(fasta):
    return
  if cram_path.ends_with(".cram"):
    stderr.write_line("[mosdepth] ERROR: specify a reference file (or set REF_PATH env var) for decoding CRAM")
    quit(1)

when(isMainModule):
  when not defined(release) and not defined(lto):
    stderr.write_line "[mosdepth] WARNING: built in debug mode; will be slow"

  let version = "mosdepth 0.2.6"
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  $version

  Usage: mosdepth [options] <prefix> <BAM-or-CRAM>

Arguments:

  <prefix>       outputs: `{prefix}.mosdepth.dist.txt`
                          `{prefix}.mosdepth.summary.txt`
                          `{prefix}.per-base.bed.gz` (unless -n/--no-per-base is specified)
                          `{prefix}.regions.bed.gz` (if --by is specified)
                          `{prefix}.quantized.bed.gz` (if --quantize is specified)
                          `{prefix}.thresholds.bed.gz` (if --thresholds is specified)

  <BAM-or-CRAM>  the alignment file for which to calculate depth.

Common Options:

  -t --threads <threads>     number of BAM decompression threads [default: 0]
  -c --chrom <chrom>         chromosome to restrict depth calculation.
  -b --by <bed|window>       optional BED file or (integer) window-sizes.
  -n --no-per-base           dont output per-base depth. skipping this output will speed execution
                             substantially. prefer quantized or thresholded values if possible.
  -f --fasta <fasta>         fasta file for use with CRAM files [default: $env_fasta].

Other options:

  -F --flag <FLAG>              exclude reads with any of the bits in FLAG set [default: 1796]
  -i --include-flag <FLAG>      only include reads with any of the bits in FLAG set. default is unset. [default: 0]
  -x --fast-mode                dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases).
  -q --quantize <segments>      write quantized output see docs for description.
  -Q --mapq <mapq>              mapping quality threshold [default: 0]
  -T --thresholds <thresholds>  for each interval in --by, write number of bases covered by at
                                least threshold bases. Specify multiple integer values separated
                                by ','.
  -m --use-median               output median of each region (in --by) instead of mean.
  -R --read-groups <string>     only calculate depth for these comma-separated read groups IDs.
  -h --help                     show help
  """ % ["version", version, "env_fasta", env_fasta])

  var args: Table[string, Value]
  try:
    args = docopt(doc, version = version, quit=false)
  except DocoptExit:
    echo (ref DocoptExit)(get_current_exception()).usage
    quit "error parsing arguments"

  let mapq = S.parse_int($args["--mapq"])
  var
    region: string
    thresholds: seq[int] = threshold_args($args["--thresholds"])
    fast_mode:bool = args["--fast-mode"]
    use_median:bool = args["--use-median"]

  if $args["--by"] != "nil":
    region = $args["--by"]
  else:
    if thresholds.len != 0:
      stderr.write_line("[mosdepth] error --thresholds can noly be used when --by is specified.")
      quit(2)
  GC_disableMarkAndSweep()
  var fasta: cstring = nil
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    eflag: uint16 = uint16(S.parse_int($args["--flag"]))
    iflag: uint16 = uint16(S.parse_int($args["--include-flag"]))
    threads = S.parse_int($args["--threads"])
    chrom = region_line_to_region($args["--chrom"])
    bam:Bam
  check_cram_has_ref($args["<BAM-or-CRAM>"], $args["--fasta"])
  open(bam, $args["<BAM-or-CRAM>"], threads=threads, index=true, fai=fasta)
  if bam.idx == nil:
    stderr.write_line("[mosdepth] error alignment file must be indexed")
    quit(2)

  var opts = SamField.SAM_FLAG.int or SamField.SAM_RNAME.int or SamField.SAM_POS.int or SamField.SAM_MAPQ.int or SamField.SAM_CIGAR.int
  if not fast_mode:
      opts = opts or SamField.SAM_QNAME.int or SamField.SAM_RNEXT.int or SamField.SAM_PNEXT.int #or SamField.SAM_TLEN.int

  if $args["--read-groups"] != "nil":
    opts = opts or SamField.SAM_RGAUX.int

  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, opts)
  discard bam.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)
  check_chrom(chrom, bam.hdr.targets)

  main(bam, chrom, mapq, eflag, iflag, region, thresholds, fast_mode, args, use_median=use_median)
