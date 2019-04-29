#!/bin/bash


#NOTE that to use this outside of travis, you'll want to export BRANCH=master

set -euo pipefail
export base=$(pwd)

BRANCH=devel

sudo apt-get update
sudo apt-get -qy install bwa make build-essential cmake libncurses-dev ncurses-dev libbz2-dev lzma-dev liblzma-dev \
     curl  libssl-dev libtool autoconf automake libcurl4-openssl-dev

git clone -b $BRANCH --depth 1 git://github.com/nim-lang/nim nim-$BRANCH/
cd nim-$BRANCH
sh build_all.sh

export PATH=$PATH:$base/nim-$BRANCH/bin/
cd $base
nimble refresh
$base/nim-$BRANCH/bin/nimble install -y

git clone --recursive https://github.com/samtools/htslib.git

cd htslib && git checkout 1.9 && autoheader && autoconf && ./configure --enable-libcurl
cd ..
make -j 4 -C htslib
export LD_LIBRARY_PATH=$base/htslib
ls -lh $base/htslib/*.so
