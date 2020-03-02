# Package

version       = "0.2.9"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.1", "docopt >= 0.6.8", "nim >= 1.0.0"

bin = @["mosdepth"]
skipDirs = @["tests"]
skipFiles = @["GT04008021.bam"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

