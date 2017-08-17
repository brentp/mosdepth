import hts as hts
import tables
import strutils as S
import os
import docopt
#import binaryheap as bh

proc dump(arr: seq[uint16], chrom: string) =
  var last_start = 0
  var last_cov = uint16(0)
  for i, cov in pairs(arr):
    if cov != last_cov:
      if last_cov == 0:
        last_cov = cov
        last_start = i
        continue
      echo chrom, "\t", intToStr(last_start), "\t", intToStr(i), "\t", $last_cov
      last_start = i
      last_cov = cov

proc main(path: string, threads:int) =
  var bam = hts.open_hts(path)
  var seqs = bam.hdr.targets

  var tgt: hts.Target
  var arr: seq[uint16]

  #var heap = bh.new_heap[hts.Range]() do (a, b: Range) -> int:
  #  return a.start - b.start
  for rec in bam:
    if rec.flag.unmapped: continue
    if tgt == nil or tgt.tid != rec.b.core.tid:
        if tgt != nil:
          dump(arr, tgt.name)
        tgt = seqs[rec.b.core.tid]
        arr = new_seq[uint16](tgt.length)
    for piece in rec.cigar.ref_coverage(rec.start):
      for i in piece.start..<piece.stop:
        arr[i]+=1
  if tgt != nil:
    dump(arr, tgt.name)


when(isMainModule):

  let doc = """
  nimdepth

  Usage: nimdepth [options] BAM
  
  -t --threads <threads>  number of threads to use [default: 1]
  -h --help               show help
  """

  let args = docopt(doc, version = "nimdepth 0.1.1", quit=true)


  main($args["BAM"], S.parse_int($args["--threads"]))
