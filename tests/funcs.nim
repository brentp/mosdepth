import unittest
import mosdepth
import depthstat
import os

suite "mosdepth-suite":

  test "depthstat min":
      var d = newSeq[int32]()
      var t = newDepthStat(d)
      check t.min_depth > 0

      var dd: depth_stat
      # "not that this always starts as 0 so we must set clear to increase it"
      check dd.min_depth == 0
      dd.clear()
      check dd.min_depth > 0

  test "test-quantize-args":

    var rs = get_quantize_args(":1")
    check rs[0] == 0
    check rs[1] == 1

    rs = get_quantize_args("0:1:4:")
    check rs[0] == 0
    check rs[1] == 1
    check rs[2] == 4
    check rs[3] == high(int)

  test "linear-search":

    var bins = @[10, 22, 44, 99]
    var idx: int

    for i, v in bins:
      linear_search(v, bins, idx.addr)
      check idx == i
    linear_search(8, bins, idx.addr)
    check idx == -1

    linear_search(800, bins, idx.addr)
    check idx == -1

    bins = @[0, 1, high(int)]
    linear_search(0, bins, idx.addr)
    check idx == 0
    linear_search(-1, bins, idx.addr)
    check idx == -1

    linear_search(99999, bins, idx.addr)
    check idx == 1

  test "lookup":
    var bins = @[10, 22, 44, 99]
    var m = make_lookup(bins)
    check m[0] == "10:22"
    check m[1] == "22:44"
    check m[2] == "44:99"
    check m.len == 3

    bins = @[0, 10]
    m = make_lookup(bins)
    check m[0] == "0:10"
    check m.len == 1

    bins = get_quantize_args("0:1:4:")
    m = make_lookup(bins)

    check m.len == 3
    check m[2] == "4:inf"

  test "threshold-args":

    var ts = threshold_args("1,2,3")
    check ts[0] == 1
    check ts[1] == 2
    check ts[2] == 3
    check ts.len == 3

  test "name splitting":

    var r = region_line_to_region("Super-Scaffold_52")
    check r.chrom == "Super-Scaffold_52"
    check r.start == 0
    check r.stop == 0

    r = region_line_to_region("Super-Scaffold_52:2-1000")
    check r.chrom == "Super-Scaffold_52"
    check r.start == 1
    check r.stop == 1000
