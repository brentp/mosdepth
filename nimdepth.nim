import hts as hts
import tables
import strutils as S
import os
import docopt
import tables

proc setvbuf(stream: File, buf: cstring, buftype: int, size: int32): int {.importc: "setvbuf", header:"<stdio.h>".}

proc dump(arr: seq[int16], chrom: string) =
  var
    last_start = -1
    last_depth = 0
    depth = 0
  
  for i, change in pairs(arr):
    depth += int(change)
    if depth == last_depth: continue
    if last_depth != 0:
      stdout.write_line chrom & "\t" & intToStr(last_start) & "\t" & intToStr(i) & "\t" & $last_depth

    last_start = i
    last_depth = depth

proc inc_coverage*(c: Cigar, ipos: int = 0, arr: var seq[int16]) {. inline .} =
  if c.len == 1 and c[0].op == CigarOp(match):
    inc(arr[ipos])
    dec(arr[ipos + c[0].len])
    return

  var pos = ipos
  var last_stop = 0
  for op in c:
    if not op.consumes_reference:
      continue
    var olen = op.len
    if op.consumes_query:
      if pos != last_stop:
        inc(arr[pos]) # increment the incoming pos
        if last_stop != 0:
          dec(arr[last_stop]) # and the last stop before overwriting
        last_stop = pos + olen
      else:
        last_stop = pos + olen
    pos += olen
  if last_stop != 0:
    dec(arr[last_stop])

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
  var bam = hts.open_hts(path, threads=threads, index=true)
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
          dec(arr[rec.start])
          inc(arr[mate.stop])
          #assert rec.start < mate.stop
          #stderr.write_line rec, " -> ", rec.cigar
          #stderr.write_line mate, " -> ", rec.cigar
          #stderr.write_line rec.start, "...", mate.stop
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
