# Package

version       = "0.3.3"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.21", "docopt >= 0.6.8", "nim >= 1.0.0", "https://github.com/brentp/d4-nim"

bin = @["mosdepth"]
skipDirs = @["tests"]
skipFiles = @["GT04008021.bam"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

