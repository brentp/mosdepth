# Package

version       = "0.1.9"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.0", "hts >= 0.1.0", "docopt"

bin = @["mosdepth"]
skipDirs = @["tests"]

task osx_build, "build binary for osx":
  exec "nim c -o:mosdepth_osx --os:macosx --cpu:amd64 --compile_only --gen_script -c mosdepth.nim"

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

