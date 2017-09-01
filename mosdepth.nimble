# Package

version       = "0.1.4"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.0", "hts >= 0.1.0", "docopt"

bin = @["mosdepth"]

task osx_build, "build binary for osx":
  exec "nim c -o:mosdepth_osx --os:macosx --cpu:amd64 --compile_only --gen_script -c mosdepth.nim"
