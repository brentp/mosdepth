import std/[os, osproc, streams, strformat, tempfiles, unittest]
import hts/bam
import d4

suite "D4 per-base output":
  test "coverage sentinel is excluded and high depths survive":
    let dir = createTempDir("mosdepth-d4-", "")
    defer: removeDir(dir)
    # Include values near dictionary boundaries and high depths to detect data loss.
    let depths = @[127'i32, 128, 254, 255, 256, 257, 5000]
    let bam = dir / "depths.bam"
    let prefix = dir / "depths"
    var header = Header()
    header.from_string("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:7\n")
    var writer: Bam
    doAssert writer.open(bam.cstring, mode = "w")
    writer.write_header(header)
    var record = NewRecord(header)
    # Write one-base alignments to produce the expected depth at each position.
    for pos, depth in depths:
      for read in 0..<depth.int:
        record.from_string(&"read-{pos}-{read}\t0\tchr1\t{pos + 1}\t60\t1M\t*\t0\t0\tA\tI")
        writer.write(record)
    writer.close()
    xam_index(bam)
    # Run mosdepth on the generated BAM to produce D4 output.
    let process = startProcess(getAppDir().parentDir / "mosdepth",
      args = ["--d4", "--fast-mode", prefix, bam], options = {poStdErrToStdOut})
    let output = process.outputStream.readAll()
    doAssert process.waitForExit() == 0, output
    process.close()
    # Read back the D4 output and verify the depth at every reference base.
    var track: D4
    doAssert track.open(prefix & ".per-base.d4")
    check track.values("chr1") == depths
    track.close()
