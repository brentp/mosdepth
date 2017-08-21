import hts as hts
import tables
import strutils as S
import algorithm as alg
import sequtils as sequtils
import os
import docopt
import tables

proc setvbuf(stream: File, buf: cstring, buftype: int, size: int32): int {.importc: "setvbuf", header:"<stdio.h>".}

proc dump(arr: seq[int16], chrom: string) =
  var
    last_start = -1
    last_depth = -1
    depth = 0
    i = -1
  
  for change in arr:
    inc(i)
    depth += int(change)
    if depth == last_depth: continue
    if last_depth != -1:
      stdout.write_line chrom & "\t" & intToStr(i) & "\t" & $last_depth

    last_start = i
    last_depth = depth
  if last_depth != -1:
    stdout.write_line chrom & "\t" & intToStr(i) & "\t" & $last_depth

type
  pair = tuple[pos: int, value: int16]

proc pair_sort(a, b: pair): int =
   return a.pos - b.pos

iterator gen_start_ends(c: Cigar, ipos: int): pair =
  if c.len == 1 and c[0].op == CigarOp(match):
    yield (ipos, int16(1))
    yield (ipos + c[0].len, int16(-1))
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
          yield (pos, int16(1))
          if last_stop != 0:
            yield (last_stop, int16(-1))
        last_stop = pos + olen
      pos += olen
    if last_stop != 0:
      yield (last_stop, int16(-1))

proc inc_coverage(c: Cigar, ipos: int = 0, arr: var seq[int16]) {. inline .} =
  for p in gen_start_ends(c, ipos):
    arr[p.pos] += p.value

iterator regions(bam: Bam, region: string): Record =
  if region == nil:
    for r in bam:
      yield r
  else:
    for r in bam.querys(region):
      yield r

proc main(path: string, threads:int=0, mapq:int= -1, region: string = "") =
  GC_disableMarkAndSweep()
  discard setvbuf(stdout, nil, 0, 16384)
  var bam = hts.open_hts(path, threads=threads, index=region != nil)
  var seqs = bam.hdr.targets

  var tgt: hts.Target
  var arr: seq[int16]
  var mate: Record
  #var heap = bh.new_heap[hts.Range]() do (a, b: Range) -> int:
  #  return a.start - b.start
  var seen = newTable[string, Record]()
  for rec in bam.regions(region):
    if rec.flag.has_flag(BAM_FUNMAP or BAM_FQCFAIL or BAM_FDUP or BAM_FSECONDARY): continue
    if int(rec.qual) < mapq: continue
    if tgt == nil or tgt.tid != rec.b.core.tid:
        if tgt != nil:
          dump(arr, tgt.name)
        tgt = seqs[rec.b.core.tid]
        arr = new_seq[int16](tgt.length)
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
              if p.value == -1 and pair_depth == 2:
                #if len(ses) > 4: stderr.write_line last_pos, " ", p.pos
                dec(arr[last_pos])
                inc(arr[p.pos])
              pair_depth += p.value
              last_pos = p.pos
            assert pair_depth == 0

    inc_coverage(rec.cigar, rec.start, arr)

  if tgt != nil:
    dump(arr, tgt.name)

when(isMainModule):

  let doc = """
  nimdepth

  Usage: nimdepth [options] <BAM>
  
  -t --threads <threads>  number of threads to use [default: 0]
  -Q --mapq <mapq>        mapping quality threshold [default: 0]
  -r --region <region>    region to restrict depth calc.
  -h --help               show help
  """

  let args = docopt(doc, version = "nimdepth 0.1.1")
  var mapq = S.parse_int($args["--mapq"])

  main($args["<BAM>"], S.parse_int($args["--threads"]), mapq, region=($args["--region"]))
