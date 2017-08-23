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
  depth_t = tuple[start: uint32, stop: uint32, value: uint32, tid: uint32]
  region_t = ref object
    chrom: string
    start: int
    stop: int

proc length(r: region_t): int =
  return r.stop - r.start

proc `$`(r: region_t): string =
  if r == nil:
    return nil
  if r.stop != 0:
    return format("$1:$2-$3", r.chrom, r.start + 1, r.stop)
  else:
    return format("$1:$2", r.chrom, r.start + 1)

iterator gen_depths(arr: seq[int32], offset: uint32, istop: int, tid: uint32): depth_t =
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
    var stop = region.stop
    if stop == 0:
      for tgt in bam.hdr.targets:
        if tgt.name == region.chrom:
          stop = int(tgt.length)
          break
    if stop == 0:
      quit(format("chromosome: $1 not found", region.chrom))
    for r in bam.query(region.chrom, region.start, stop):
      yield r

proc bed_line_to_region(line: string): region_t =
   var cse = sequtils.to_seq(line.strip().split("\t"))
   var s = S.parse_int(cse[1])
   var e = S.parse_int(cse[2])
   var reg = region_t(chrom: cse[0], start:s, stop: e)
   return reg

proc region_line_to_region(region: string): region_t =
  if region == nil or region == "" or region == "nil":
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

iterator depth(bam: hts.Bam, arr: var seq[int32], region: var region_t, mapq:int= -1): depth_t =
  var seqs = bam.hdr.targets

  var tgt: hts.Target
  var mate: Record
  var seen = newTable[string, Record]()

  for rec in bam.regions(region):
    if int(rec.qual) < mapq: continue
    if rec.flag.has_flag(BAM_FUNMAP or BAM_FQCFAIL or BAM_FDUP or BAM_FSECONDARY): continue
    if tgt == nil or tgt.tid != rec.b.core.tid:
        if tgt != nil:
          for p in gen_depths(arr, 0, 0, uint32(tgt.tid)): yield p
          flushFile(stdout)
        tgt = seqs[rec.b.core.tid]
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
    for p in gen_depths(arr, 0, 0, uint32(tgt.tid)): yield p
  flushFile(stdout)

proc bed_main(path: string, bed: string, threads:int=0, mapq:int= -1) =
  var bam = hts.open_hts(path, threads=1, index=true)
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  var arr: seq[int32]
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    var rl = bed_line_to_region($kstr.s)
    #main(bam, arr, rl, mapq)

  hts.free(kstr.s)

when(isMainModule):

  let doc = """
  nimdepth

  Usage: nimdepth [options] <BAM>
  
  -t --threads <threads>  number of threads to use [default: 0]
  -c --chrom <chrom>      chromosome to restrict depth calculation.
  -Q --mapq <mapq>        mapping quality threshold [default: 0]
  -L --bed <bed>          BED file of regions for which to calculate depth.
  -h --help               show help
  """

  let args = docopt(doc, version = "nimdepth 0.1.1")
  let mapq = S.parse_int($args["--mapq"])
  var bed_based = false 
  if $args["--bed"] != "nil":
    bed_based = true
  GC_disableMarkAndSweep()
  discard setvbuf(stdout, nil, 0, 16384)

  var
    arr:seq[int32]
    threads = S.parse_int($args["--threads"])
    chrom = region_line_to_region($args["--chrom"])
    bam = hts.open_hts($args["<BAM>"], threads=threads, index=chrom != nil)
    targets = bam.hdr.targets()
    last_tid = uint32(0)
    target = targets[int(last_tid)].name & "\t"

  for region in depth(bam, arr, chrom, mapq):
    if region.tid != last_tid:
      last_tid = region.tid
      target = targets[int(last_tid)].name & "\t"
    #stdout.write_line(target & intToStr(int(region.start)) & "\t" & intToStr(int(region.stop)) & "\t" & intToStr(int(region.value)))
    #stdout.write_line(target & intToStr(int(region.stop)) & "\t" & intToStr(int(region.value)))
    stdout.write_line(intToStr(int(region.stop)) & "\t" & intToStr(int(region.value)))
