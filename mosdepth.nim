import hts as hts
import tables
import strutils as S
import algorithm as alg
import sequtils as sequtils
import strutils as su
import os
import docopt
import tables
import times

proc setvbuf(stream: File, buf: cstring, buftype: int, size: int32): int {.importc: "setvbuf", header:"<stdio.h>".}
type
  pair = tuple[pos: int, value: int32]
  depth_t = tuple[start: uint32, stop: uint32, value: uint32, tid: uint32]
  region_t = ref object
    chrom: string
    start: uint32
    stop: uint32
    name: string

  coverage_t = seq[int32]

proc `$`(r: region_t): string =
  if r == nil:
    return nil
  if r.stop != 0:
    return format("$1:$2-$3", r.chrom, r.start + 1, r.stop)
  else:
    return format("$1:$2", r.chrom, r.start + 1)

proc to_coverage(c: var coverage_t) =
  # to_coverage converts from an array of start/end inc/decs to actual coverage.
  var d = int32(0)
  for i, v in pairs(c):
    d += v
    c[i] = d

proc length(r: region_t): int =
  return int(r.stop - r.start)

iterator gen_depths(arr: coverage_t, offset: uint32, istop: int, tid: uint32): depth_t =
  # given `arr` with values in each index indicating the number of reads
  # starting or ending at that location, generate depths.
  # offset is only used for a region like chr6:200-30000, in which case, offset will be 200
  var
    last_depth = -1
    depth = 0
    i = uint32(0)
    last_i = uint32(0)
    stop: uint32
    last_yield: bool
  if istop <= 0:
    stop = uint32(len(arr)-1)
  else:
    stop = uint32(istop)
  # even with an offset, have to start from the beginning of the array
  # to get the proper depth.
  for change in arr:
    depth += int(change)
    if i == stop:
      break
    last_yield = false
    if i < offset or depth == last_depth:
      inc(i)
      continue

    if last_depth != -1:
      yield (last_i, i, uint32(last_depth), tid)
      last_yield = true

    last_depth = depth
    last_i = i
    inc(i)
  if last_i < uint32(len(arr)-1):
    yield (last_i, uint32(len(arr)-1), uint32(0), tid)

  # this is horrible, but it works. we don't know
  # if we've already printed the record on not.
  if  last_yield != (last_i == i):
    if last_i == i:
      yield (last_i - 1, i, uint32(last_depth), tid)
    else:
      yield (last_i, i, uint32(last_depth), tid)

proc pair_sort(a, b: pair): int =
   return a.pos - b.pos

iterator gen_start_ends(c: Cigar, ipos: int): pair =
  # generate start, end pairs given a cigar string and a position offset.
  if c.len == 1 and c[0].op == CigarOp(match):
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

proc inc_coverage(c: Cigar, ipos: int = 0, arr: var seq[int32]) {. inline .} =
  for p in gen_start_ends(c, ipos):
      arr[p.pos] += p.value

iterator regions(bam: hts.Bam, region: region_t, tid: int, targets: seq[hts.Target]): Record =
  if region == nil:
    for r in bam:
      yield r
  elif region != nil:
    var stop = region.stop
    if tid != -1:
      if stop == 0:
        stop = targets[tid].length
      for r in bam.queryi(uint32(tid), int(region.start), int(stop)):
        yield r
    else:
      for r in bam.query(region.chrom, int(region.start), int(stop)):
        yield r

proc bed_line_to_region(line: string): region_t =
   var
     cse = sequtils.to_seq(line.strip().split("\t"))

   if len(cse) < 3:
     stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
     return nil

   var
     s = S.parse_int(cse[1])
     e = S.parse_int(cse[2])
     reg = region_t(chrom: cse[0], start: uint32(s), stop: uint32(e))
   if len(cse) > 3:
     reg.name = cse[3]
   return reg

proc region_line_to_region(region: string): region_t =
  if region == nil or region == "" or region == "nil":
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

iterator coverage(bam: hts.Bam, arr: var coverage_t, region: var region_t, mapq:int= -1, eflag: uint16=1796): int =
  # depth updates arr in-place and yields the tid for each chrom.
  var
    targets = bam.hdr.targets
    tgt: hts.Target
    mate: Record
    seen = newTable[string, Record]()

  var tid = if region != nil: get_tid(targets, region.chrom) else: -1

  for rec in bam.regions(region, tid, targets):
    if int(rec.qual) < mapq: continue
    if (rec.flag and eflag) != 0:
      continue
    if tgt == nil or tgt.tid != rec.b.core.tid:
        if tgt != nil:
          yield tgt.tid
          #for p in gen_depths(arr, 0, 0, uint32(tgt.tid)): yield p
          flushFile(stdout)
        tgt = targets[rec.b.core.tid]
        if arr == nil or len(arr) != int(tgt.length+1):
          # must create a new array in some cases.
          if arr == nil or len(arr) < int(tgt.length+1):
            arr = new_seq[int32](tgt.length+1)
          else:
            # otherwise can re-use and zero
            arr.set_len(int(tgt.length+1))
            for i in 0..<len(arr):
              arr[i] = 0

        seen.clear()
    # rec:   --------------
    # mate:             ------------
    # handle overlapping mate pairs.
    if rec.flag.proper_pair:
      if rec.b.core.tid == rec.b.core.mtid and rec.stop > rec.matepos and rec.start < rec.matepos:
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
            assert pair_depth == 0, $rec.qname & ":" & $rec & " " & $mate.qname & ":" & $mate & " " & $pair_depth
    inc_coverage(rec.cigar, rec.start, arr)

  if tgt != nil:
    yield tgt.tid
    #for p in gen_depths(arr, 0, 0, uint32(tgt.tid)): yield p
  flushFile(stdout)

iterator bed_gen(bed: string): region_t =
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    yield bed_line_to_region($kstr.s)

  hts.free(kstr.s)

iterator window_gen(window: uint32, targets: seq[hts.Target]): region_t =
  for t in targets:
    var start:uint32 = 0
    while start + window < t.length:
      yield region_t(chrom: t.name, start: start, stop: start + window)
      start += window
    if start != t.length:
      yield region_t(chrom: t.name, start: start, stop: t.length)

iterator region_gen(bed_or_window: string, targets: seq[hts.Target]): region_t =
  if bed_or_window.isdigit():
    for r in window_gen(uint32(S.parse_int(bed_or_window)), targets): yield r
  else:
    for r in bed_gen(bed_or_window): yield r

proc imean(vals: coverage_t, start:uint32, stop:uint32): float64 =
  if vals == nil or start > uint32(len(vals)):
    return 0
  for i in start .. <stop:
    if int(i) == len(vals): break
    result += float64(vals[int(i)])
  result /= float64(stop-start)

const MAX_COVERAGE = int32(400000)

proc inc(d: var seq[int32], coverage: coverage_t, start:uint32, stop:uint32) =
  var v:int32
  var L = int32(len(d)-1)
  for i in start .. <stop:
    v = coverage[int(i)]
    if v > MAX_COVERAGE:
      v = MAX_COVERAGE - 10
    if v >= L:
      d.set_len(v + 10)
      L = int32(len(d)-1)
    if v < 0: continue
    inc(d[v])

proc write_distribution(d: var seq[int32], path:string) =
  var fh:File
  if not open(fh, path, fmWrite):
    stderr.write_line("[mosdepth] could not open file:", path)
    quit(1)
  var sum: int64
  for v in d: sum += int64(v)
  var cum: float64 = 0
  # reverse and then cumsum so that e.g. a value of 1 is the proportion of
  # bases with a coverage of at least 1.
  reverse(d)
  # skip until we see a non-zero value for high-coverage end of array.
  for i, v in pairs(d):
    var irev = len(d) - i - 1
    if irev > 300 and v == 0: continue
    cum += float64(v) / float64(sum)
    if cum < 5e-6: continue
    fh.write_line($irev & "\t" & su.format_float(cum, ffDecimal, precision=5))
  fh.close()

proc get_targets(targets: seq[hts.Target], r: region_t): seq[hts.Target] =
  if r == nil:
    return targets
  result = new_seq[hts.Target](1)
  for t in targets:
    if t.name == r.chrom:
      result[0] = t
      return result

proc window_main(bam: hts.Bam, chrom: region_t, mapq: int, eflag: uint16, args: Table[string, docopt.Value]) =
  var targets = bam.hdr.targets
  # windows are either from regions, or fixed-length windows.
  # we assume the input is sorted by chrom.
  var
    sub_targets = get_targets(targets, chrom)
    distribution: seq[int32]
    last_chrom = ""
    rchrom : region_t
    found = false
    target: string
    arr: coverage_t

  if $args["--distribution"] != "nil":
    distribution = new_seq[int32](1000)

  for r in region_gen($args["--by"], sub_targets):
    if r == nil: continue
    if chrom != nil and r.chrom != chrom.chrom: continue
    # the firs time seeing the chrom, we fill the coverage array for
    # the entire chrom.
    if r.chrom != last_chrom:
      target = r.chrom & "\t"
      var j = 0
      rchrom = region_t(chrom: r.chrom)
      for tid in coverage(bam, arr, rchrom, mapq, eflag):
        arr.to_coverage()
        inc(j)
      last_chrom = r.chrom
      if j == 0: # didn't find this chrom
        stderr.write_line "[mosdepth] chromosome: ", r.chrom, " not found in alignments"
        found = false
      else:
        found = true
    #if not found: continue # now just writes from imean

    var me = imean(arr, r.start, r.stop)
    var m = su.format_float(me, ffDecimal, precision=2)
    if r.name == nil:
      stdout.write_line(target, intToStr(int(r.start)), "\t", intToStr(int(r.stop)), "\t", m)
    else:
      stdout.write_line(target, intToStr(int(r.start)), "\t", intToStr(int(r.stop)), "\t", r.name, "\t", m)
    if distribution != nil and arr != nil and found:
      distribution.inc(arr, r.start, r.stop)
  if distribution != nil:
    write_distribution(distribution, $args["--distribution"])

proc check_chrom(r: region_t, targets: seq[Target]) =
  if r == nil: return
  for t in targets:
    if t.name == r.chrom:
      return
  stderr.write_line "[mosdepth] chromosome ", r.chrom, " not found"
  quit(1)

when(isMainModule):
  when not defined(release):
    stderr.write_line "[mosdepth] WARNING: built debug mode. will be slow"

  let doc = """
  mosdepth

  Usage: mosdepth [options] <BAM-or-CRAM>

Common Options:
  
  -t --threads <threads>     number of BAM decompression threads [default: 0]
  -c --chrom <chrom>         chromosome to restrict depth calculation.
  -b --by <bed|window>       BED file of regions or an (integer) window-size.
  -d --distribution <file>   a cumulative distribution file (coverage, proportion).
  -f --fasta <fasta>         fasta file for use with CRAM files.

Other options:

  -F --flag <FLAG>           exclude reads with any of the bits in FLAG set [default: 1796]
  -Q --mapq <mapq>           mapping quality threshold [default: 0]
  -h --help                  show help
  """

  let args = docopt(doc, version = "mosdepth 0.1.7")
  let mapq = S.parse_int($args["--mapq"])
  var window_based = false 
  if $args["--by"] != "nil":
    window_based = true
  GC_disableMarkAndSweep()
  discard setvbuf(stdout, nil, 0, 16384)
  var fasta: cstring = nil
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    arr:coverage_t
    eflag: uint16 = uint16(S.parse_int($args["--flag"]))
    threads = S.parse_int($args["--threads"])
    chrom = region_line_to_region($args["--chrom"])
    bam = hts.open_hts($args["<BAM-or-CRAM>"], threads=threads, index=chrom != nil or window_based, fai=fasta)
    targets = bam.hdr.targets()
    last_tid = uint32(0)
    target = targets[int(last_tid)].name & "\t"

  discard bam.set_fields(SamField.SAM_QNAME, SamField.SAM_FLAG, SamField.SAM_RNAME,
                         SamField.SAM_POS, SamField.SAM_MAPQ, SamField.SAM_CIGAR,
                         SamField.SAM_RNEXT, SamField.SAM_PNEXT, SamField.SAM_TLEN,
                         SamField.SAM_QUAL, SamField.SAM_AUX)
  discard bam.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)

  check_chrom(chrom, targets)

  if not window_based:
    var distribution: seq[int32]
    if $args["--distribution"] != "nil":
      distribution = new_seq[int32](1000)
    for tid in coverage(bam, arr, chrom, mapq, eflag):
        target = targets[int(tid)].name & "\t"
        for p in gen_depths(arr, 0, 0, uint32(tid)):
            stdout.write_line(target, intToStr(int(p.stop)), "\t", intToStr(int(p.value)))
        if distribution != nil:
          arr.to_coverage()
          distribution.inc(arr, uint32(0), uint32(len(arr)))
    if distribution != nil:
      write_distribution(distribution, $args["--distribution"])
    quit()

  window_main(bam, chrom, mapq, eflag, args)
