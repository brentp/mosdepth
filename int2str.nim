import bitops

#[

this intToStr is adpated from: https://github.com/miloyip/itoa-benchmark countlut method.

https://github.com/miloyip/itoa-benchmark/blob/master/license.txt
Copyright (C) 2014 Milo Yip

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
]#


const gDigitsLut = [
    '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
    '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
    '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
    '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
    '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
    '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
    '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
    '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
    '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
    '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9'
];

const powers_of_10 = [
        0'i32,
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
        10000000,
        100000000,
        1000000000]

proc countdigits(value:int32): int {.inline.} =
  result = (32 - countLeadingZeroBits(value or 1)) * 1233 shr 12
  result = result - int(value < powers_of_10[result]) + 1

proc fastIntToStr*(value:int32, outstr:var string, offset:int=0) {.inline.} =
  outstr.setLen(offset + countdigits(value))
  var value = value
  var L = outstr.high
  while value >= 100:
    let i = (value mod 100) shl 1
    value = int32(value / 100)

    outstr[L] = gDigitsLut[i + 1]
    outstr[L-1] = gDigitsLut[i]
    L -= 2
  if value < 10:
    outstr[L] = (value + '0'.int).char

  else:
    let i = value shl 1
    outstr[L] = gDigitsLut[i + 1]
    outstr[L-1] = gDigitsLut[i]


when isMainModule:

  import times
  import strutils

  var t0 = cpuTime()
  var outstr = newString(10)
  for i in 0'i32..200_000_000:
    fastIntToStr(i, outstr)
    doAssert outstr == $i
  echo cpuTime() - t0

  t0 = cpuTime()
  echo "only to 100m for $"
  for i in 0'i32..100_000_000:
    let outstr = intToStr(i)
  echo cpuTime() - t0
