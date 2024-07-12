# Package

version       = "0.3.8"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.22", "docopt == 0.7.1", "nim >= 1.0.0", "https://github.com/brentp/d4-nim >= 0.0.3"

bin = @["mosdepth"]
skipDirs = @["tests"]
skipFiles = @["GT04008021.bam"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

