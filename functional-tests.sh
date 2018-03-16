#!/bin/bash
unset MOSDEPTH_Q0
unset MOSDEPTH_Q1
unset MOSDEPTH_Q2
unset MOSDEPTH_PRECISION

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

set -e
nim c --boundChecks:off mosdepth.nim
set +e
exe=./mosdepth
bam=/data/human/NA12878.subset.bam


unset REF_PATH
run cram_no_ref $exe xx tests/tt.cram
assert_exit_code 1

run overlapM $exe t tests/ovl.bam
assert_exit_code 0
assert_equal "$(zgrep ^MT t.per-base.bed.gz)" "MT	0	80	1
MT	80	16569	0"
assert_equal "$(zgrep -w ^1 t.per-base.bed.gz)" "1	0	249250621	0"

run missing_chrom $exe -c nonexistent --by 20000 t tests/ovl.bam
assert_in_stderr "[mosdepth] chromosome nonexistent not found"
assert_exit_code 1

run unordered_bed $exe --by tests/unordered.bed t tests/ovl.bam
assert_exit_code 0
assert_equal $(zcat t.regions.bed.gz | wc -l) 2

# theres data left in the bam but the region tree is empty...
run missing_bed_chrom $exe --by tests/missing.bed t tests/ovl.bam
assert_exit_code 0

run big_window $exe t tests/ovl.bam --by 100000000
assert_exit_code 0
assert_equal $(zgrep -c "MT" t.per-base.bed.gz) 2
assert_equal "MT	0	16569	0.00" "$(zgrep ^MT t.regions.bed.gz)"

unset MOSDEPTH_Q0
unset MOSDEPTH_Q1
unset MOSDEPTH_Q2
rm -f t.quantized.bed.gz
run quantest $exe -q 0:1:1000 t tests/ovl.bam
assert_exit_code 0
assert_equal "$(zgrep ^MT t.quantized.bed.gz)" "MT	0	80	1:1000
MT	80	16569	0:1"
assert_equal "$(zgrep -w ^1 t.quantized.bed.gz)" "1	0	249250621	0:1"


rm -f t.thresholds.bed.gz*
run threshold_test $exe --by 100 -T 0,1,2,3,4,5 -c MT t tests/ovl.bam
assert_equal "$(zcat t.thresholds.bed.gz | tail -n +2 | head -1)" "MT	0	100	unknown	100	80	0	0	0	0"
assert_equal "0" "$(zcat t.thresholds.bed.gz | tail -n+2 | cut -f 7 | uniq)"
assert_exit_code 0

rm -f t.thresholds.bed.gz*
run threshold_test_by $exe --by tests/track.bed -T 0,1,2 -c MT t tests/ovl.bam
assert_equal "$(zcat t.thresholds.bed.gz | tail -n +2)" "MT	2	80	aregion	78	78	0"
assert_exit_code 0

export MOSDEPTH_Q0=AAA
export MOSDEPTH_Q1=BBB
rm -f t.quantized.bed.gz t.mosdepth.global.dist.txt
MOSDEPTH_PRECISION=7 run quantest-named-and-precision $exe -q 0:1:1000 t tests/ovl.bam
assert_exit_code 0
assert_equal "$(zgrep -w ^MT t.quantized.bed.gz)" "MT	0	80	BBB
MT	80	16569	AAA"
assert_equal "$(head -1 t.mosdepth.global.dist.txt)" "MT	1	0.0048280"

run track_header $exe --by tests/track.bed t tests/ovl.bam
assert_exit_code 0
assert_equal "$(zcat t.regions.bed.gz)" "MT	2	80	aregion	1.00"

run track_header_by $exe --by tests/bad.bed t tests/ovl.bam
assert_exit_code 1
assert_in_stderr "skipping bad bed line:MT	2"
assert_in_stderr "invalid integer: asdf"


test -e $bam || exit

run short $exe -c chrM t $bam
assert_equal "$(zgrep -w "chrM	9999	10000" t.per-base.bed.gz)" "chrM	9999	10000	19073"
assert_exit_code 0

run flag $exe -c chrM -F 4 --by 20000 tx /data/human/NA12878.subset.bam
assert_exit_code 0

