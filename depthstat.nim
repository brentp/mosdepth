
type
  depth_stat* = object
    cum_length*: int
    cum_depth*: uint64
    min_depth*: uint32
    max_depth*: uint32

proc sum64*(x: seq[int32]): uint64 =
  # Required for summing up cumulative depth
  for i in items(x): result = result + i.uint64

proc newDepthStat*(d: seq[SomeNumber]): depth_stat =
  return depth_stat(cum_length: len(d),
                    cum_depth: sum64(d),
                    min_depth: min(d).uint32,
                    max_depth: max(d).uint32)

proc clear*(ds: var depth_stat) =
    ds.cum_length = 0
    ds.cum_depth = 0
    ds.min_depth = uint32.high
    ds.max_depth = 0

proc `+`*(a, b: depth_stat): depth_stat =
    return depth_stat(cum_length: a.cum_length + b.cum_length,
                      cum_depth: a.cum_depth + b.cum_depth,
                      min_depth: min([a.min_depth, b.min_depth]),
                      max_depth: max([a.max_depth, b.max_depth]))
