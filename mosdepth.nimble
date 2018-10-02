# Package

version       = "0.2.2"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "hts >= 0.1.8", "docopt#0abba63"

bin = @["mosdepth"]
skipDirs = @["tests"]
skipFiles = @["GT04008021.bam"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

