# Package

version       = "0.2.5"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "hts >= 0.2.7", "docopt >= 0.6.8"

bin = @["mosdepth"]
skipDirs = @["tests"]
skipFiles = @["GT04008021.bam"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

