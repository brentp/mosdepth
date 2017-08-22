#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

nim c nimdepth.nim
exe=./nimdepth
bam=/data/human/NA12878.subset.bam

run short $exe -c chrM $bam
head $STDOUT_FILE
tail $STDOUT_FILE
assert_in_stdout "chrM	10000	19067"
assert_exit_code 0
