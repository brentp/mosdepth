# GitHub Copilot Instructions for mosdepth

## Project Overview
mosdepth is a fast BAM/CRAM depth calculation tool for WGS, exome, or targeted sequencing written in Nim. It provides per-base, per-window, and per-region depth calculations about 2x faster than samtools depth.

## Technologies
- **Language**: Nim (requires Nim >= 1.0.0)
- **Key Dependencies**: 
  - hts-nim (wrapper around htslib for BAM/CRAM parsing)
  - htslib (requires version 1.4+)
  - d4-nim (for D4 format support)
  - docopt (for command-line parsing)
- **Build Requirements**: 
  - Must be built with `--mm:refc` flag for optimal performance and correct operation
  - Release builds use: `nim c --mm:refc -d:danger -d:release`

## Architecture
- **Core Algorithm**: Uses a cumulative array approach where starts/stops are tracked in an array the length of each chromosome
  - Increments for read starts, decrements for read stops
  - Cumulative sum provides depth at any position
  - This is faster than pileup-based approaches but uses more memory (~32-bits * longest chromosome length)
- **Key Files**:
  - `mosdepth.nim`: Main entry point and core depth calculation logic
  - `depthstat.nim`: Statistics and distribution calculations
  - `int2str.nim`: Integer to string conversion utilities

## Coding Standards
- Follow Nim naming conventions:
  - `snake_case` for variables, functions, and types
  - Use type suffixes: `_t` for tuple types, `_s` for string tuple variants
  - Prefix for objects: use `*` for public fields
- Prefer explicit types for clarity in performance-critical code
- Use inline pragmas (`{.inline.}`) for hot-path functions
- Use shallow copies (`{.shallow.}`) for large sequences when appropriate

## Testing
- **Unit Tests**: Located in `tests/` directory
  - Run with: `nim c -r tests/all.nim`
  - Uses the standard `unittest` module
  - Test file: `tests/funcs.nim`
- **Functional Tests**: Shell-based tests in `functional-tests.sh`
  - Uses ssshtest framework
  - Tests real BAM/CRAM files with expected outputs
  - Run with: `bash ./functional-tests.sh`
- **Test Data**: Test BAM/CRAM files are in the `tests/` directory

## Building
```bash
# Development build
nimble build -Y mosdepth.nimble

# Release build
nim c --mm:refc -d:danger -d:release mosdepth.nim
```

## Key Features to Preserve
- **Mate-pair overlap detection**: Avoids double-counting overlapping mate-pairs (unless --fast-mode is used)
- **CIGAR-aware counting**: Tracks every aligned part using CIGAR operations
- **Environment variables**:
  - `MOSDEPTH_PRECISION`: Controls decimal precision in output (default: 2)
  - `MOSDEPTH_Q0`, `MOSDEPTH_Q1`, etc.: Custom labels for quantize bins
- **Output formats**: per-base BED, region BED, quantized BED, threshold BED, summary TXT, distribution TXT, D4 format

## Performance Considerations
- Memory usage is ~32-bits per base for longest chromosome (e.g., 1GB for 249MB chr1)
- Threading applies to BAM decompression only (diminishing returns after ~4 threads)
- Per-base output is expensive; prefer quantized or window-based when possible
- `--fast-mode` skips mate overlap correction and CIGAR operations for speed

## Common Pitfalls
- Must handle empty chromosomes correctly (see tests/empty-tids.bam)
- Region files must be sorted by chromosome order matching BAM header
- CRAM files require reference FASTA with `-f/--fasta` flag
- Need to handle both BAM and CRAM formats appropriately

## External Resources
- [hts-nim documentation](https://github.com/brentp/hts-nim/)
- [Nim language documentation](https://nim-lang.org/docs/)
- [D4 format specification](https://github.com/38/d4-format)
