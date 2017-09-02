#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

nim c mosdepth.nim
exe=./mosdepth
bam=/data/human/NA12878.subset.bam


run overlapM $exe tests/ovl.bam
assert_exit_code 0
assert_in_stdout "MT	80	1
MT	16569	0"


test -e $bam || exit

run short $exe -c chrM $bam
head $STDOUT_FILE
tail $STDOUT_FILE
assert_in_stdout "chrM	10000	19073"
assert_exit_code 0


run missing_chrom $exe -c nonexistent --by 20000 /data/human/NA12878.subset.bam
assert_in_stderr "[mosdepth] chromosome nonexistent not found"

assert_exit_code 1

run flag $exe -c chrM -F 4 --by 20000 /data/human/NA12878.subset.bam
assert_exit_code 0

