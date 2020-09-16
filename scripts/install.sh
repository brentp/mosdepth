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

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -qy

export HTSLIB=dynamic
sudo apt-get update
sudo apt-get install git llvm curl wget libcurl4-openssl-dev
wget https://github.com/samtools/htslib/archive/1.10.2.tar.gz
tar xzf 1.10.2.tar.gz
cd htslib-1.10.2/
autoheader && autoconf && ./configure --enable-libcurl
sudo make -j4 install
git clone https://github.com/38/d4-format
cd d4-format
~/.cargo/bin/cargo build --release
sudo cp ../d4-format/target/release/libd4binding.* /usr/local/lib
sudo cp ./d4binding/include/d4.h /usr/local/include/
sudo ldconfig


cd ..
export LD_LIBRARY_PATH=$base/htslib-1.10.2
ls -lh $base/htslib-1.10.2/*.so

nimble install -y https://github.com/brentp/d4-nim/
