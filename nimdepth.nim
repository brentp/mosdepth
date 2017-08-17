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
    if cov == last_cov: continue
    if last_cov != 0:
      echo chrom, "\t", intToStr(last_start), "\t", intToStr(i), "\t", $last_cov

    last_start = i
    last_cov = cov

proc inc_coverage*(c: Cigar, ipos: int = 0, arr: var seq[uint16], qname: string="") {. inline .} =
  if c.len == 1 and c[0].op == CigarOp(match):
    for i in ipos..<(ipos + c[0].len):
       arr[i] += 1
    return

  var pos = ipos
  var last: Range = (0, 0)
  for op in c:
    if not op.consumes_reference:
      continue
    var olen = op.len
    if op.consumes_query:
      if pos != last.stop:
        for i in last.start..<last.stop:
          arr[i] += 1
        last = (pos, pos+olen)
      else:
        last.stop = pos + olen
    pos += olen
  for i in last.start..<last.stop:
    arr[i] += 1

iterator regions(bam: Bam, region: string): Record =
  if region == nil:
    for r in bam:
      yield r
  else:
    for r in bam.querys(region):
      yield r

proc main(path: string, threads:int=0, mapq:int= -1, region: string = "") =
  var bam = hts.open_hts(path, threads=threads, index=true)
  var seqs = bam.hdr.targets

  var tgt: hts.Target
  var arr: seq[uint16]
  #var heap = bh.new_heap[hts.Range]() do (a, b: Range) -> int:
  #  return a.start - b.start
  for rec in bam.regions(region):
    if rec.flag.has_flag(BAM_FUNMAP or BAM_FQCFAIL or BAM_FDUP or BAM_FSECONDARY): continue
    if int(rec.qual) < mapq: continue
    if tgt == nil or tgt.tid != rec.b.core.tid:
        if tgt != nil:
          dump(arr, tgt.name)
        tgt = seqs[rec.b.core.tid]
        arr = new_seq[uint16](tgt.length)
    inc_coverage(rec.cigar, rec.start, arr, rec.qname)

  if tgt != nil:
    dump(arr, tgt.name)

when(isMainModule):

  let doc = """
  nimdepth

  Usage: nimdepth [options] <BAM>
  
  -t --threads <threads>  number of threads to use [default: 1]
  -Q --mapq <mapq>        mapping quality threshold [default: 0]
  -r --region <region>    region to restrict depth calc.
  -h --help               show help
  """

  let args = docopt(doc, version = "nimdepth 0.1.1")
  var mapq = S.parse_int($args["--mapq"])

  main($args["<BAM>"], S.parse_int($args["--threads"]), mapq, region=($args["--region"]))
