# Package

version       = "0.2.0"
author        = "Brent Pedersen"
description   = "fast depth"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.0", "hts >= 0.1.5", "docopt"

bin = @["mosdepth"]
skipDirs = @["tests"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

