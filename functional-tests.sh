#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

nim c mosdepth.nim
exe=./mosdepth
bam=/data/human/NA12878.subset.bam


run overlapM $exe t tests/ovl.bam
assert_exit_code 0
assert_equal "$(cat t.per-base.txt)" "MT	80	1
MT	16569	0"

run missing_chrom $exe -c nonexistent --by 20000 t tests/ovl.bam
assert_in_stderr "[mosdepth] chromosome nonexistent not found"
assert_exit_code 1

run big_window $exe t tests/ovl.bam --by 100000000
assert_exit_code 0
assert_equal $(grep -c "MT" t.per-base.txt) 2
assert_equal "MT	0	16569	0.00" "$(cat t.regions.bed)"


run track_header $exe --by tests/track.bed t tests/ovl.bam
assert_exit_code 0
assert_equal "$(cat t.regions.bed)" "MT	2	80	1.00"

run track_header $exe --by tests/bad.bed t tests/ovl.bam
assert_exit_code 1
assert_in_stderr "skipping bad bed line:MT	2"
assert_in_stderr "invalid integer: asdf"




test -e $bam || exit

run short $exe -c chrM t $bam
assert_equal "$(grep -w "chrM	10000" t.per-base.txt)" "chrM	10000	19073"
assert_exit_code 0

run flag $exe -c chrM -F 4 --by 20000 t /data/human/NA12878.subset.bam
assert_exit_code 0


