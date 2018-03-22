#!/bin/bash


#NOTE that to use this outside of travis, you'll want to export BRANCH=master

echo $(pwd)

BRANCH=${BRANCH:-master}

sudo apt-get -qy install bwa make build-essential cmake libncurses-dev ncurses-dev libbz2-dev lzma-dev liblzma-dev \
     curl  libssl-dev libtool autoconf automake libcurl4-openssl-dev

export base=$(pwd)
if [ ! -x nim-$BRANCH/bin/nim ]; then
  git clone -b $BRANCH --depth 1 git://github.com/nim-lang/nim nim-$BRANCH/
  cd nim-$BRANCH
  git clone -b $BRANCH --depth 1 git://github.com/nim-lang/csources csources/
  cd csources
  sh build.sh
  cd ..
  rm -rf csources
  bin/nim c koch
  ./koch boot -d:release
else
  cd nim-$BRANCH
  git fetch origin
  if ! git merge FETCH_HEAD | grep "Already up-to-date"; then
    bin/nim c koch
    ./koch boot -d:release
  fi
fi

export PATH=$PATH:$base/nim-$BRANCH/bin/:$PATH:$base/nimble/src
cd $base
echo $PATH


git clone --depth 1 https://github.com/nim-lang/nimble.git
cd nimble
nim c src/nimble
src/nimble install -y

set -x
cd $base
nimble refresh

echo $(which nimble)
echo $(pwd)


if [ ! -x hts-nim ]; then
    cd $base
    git clone --depth 1 https://github.com/brentp/hts-nim/
    cd hts-nim && nimble install -y
fi

set -x
cd $base
nimble install -y
git clone --recursive https://github.com/samtools/htslib.git

cd htslib && git checkout 1.5 && autoheader && autoconf && ./configure --enable-libcurl
cd ..
make -j 4 -C htslib
export LD_LIBRARY_PATH=$base/htslib
ls -lh $base/htslib/*.so
echo $base
