
type
  depth_stat* = object
    cum_length*: int
    cum_depth*: uint64
    min_depth*: uint32
    max_depth*: uint32

type CountStat*[T:SomeOrdinal] = object
  counts: seq[T]
  n: int

proc initCountStat*[T](size:int=32768): CountStat[T] =
  return CountStat[T](counts: newSeq[T](size))

proc add*[T](c: var CountStat, value: T) {.inline.} =
  c.n.inc
  if value.int > c.counts.high.int:
    c.counts[c.counts.high].inc
  elif value < 0:
    raise newException(IndexError, "error setting negative depth value:" & $value)
  else:
    c.counts[value].inc

proc median*[T](c: CountStat[T]): int {.inline.} =
  var stop_n = int(0.5 + c.n.float64 * 0.5)
  var cum = 0
  for i, cnt in c.counts:
    cum += cnt.int
    if cum >= stop_n:
      return i
  return -1

proc clear*[T](c: var CountStat[T]) {.inline.} =
  if c.n == 0: return
  c.n = 0
  zeroMem(c.counts[0].addr, sizeof(c.counts[0]) * c.counts.len)

template len*[T](c:CountStat[T]): int = c.counts.len

proc newDepthStat*[T: SomeNumber](d: seq[T]): depth_stat =
  result.cum_length = len(d)
  result.min_depth = uint32.high
  for dp in d:
    result.cum_depth += dp.uint64
    result.min_depth = min(result.min_depth, dp.uint32)
    result.max_depth = max(result.max_depth, dp.uint32)

proc clear*(ds: var depth_stat) =
    ds.cum_length = 0
    ds.cum_depth = 0
    ds.min_depth = uint32.high
    ds.max_depth = 0

proc `+`*(a, b: depth_stat): depth_stat {.inline.} =
  result = depth_stat(cum_length: a.cum_length + b.cum_length,
                      cum_depth: a.cum_depth + b.cum_depth,
                      min_depth: min(a.min_depth, b.min_depth),
                      max_depth: max(a.max_depth, b.max_depth))

when isMainModule:
  import unittest
  suite "count-stat":
    test "count-test":

      var c = initCountStat[uint32]()
      for i in 0..10:
        c.add(i.uint32)

      check c.median == 5
      check c.n == 11
      check c.len == 32768

      c.add(6)
      c.add(6)
      c.add(6)
      c.add(6)
      check c.median == 6
      c.clear()
      check c.median == 0
      check c.n == 0
