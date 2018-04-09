set -euo pipefail

NIM_VERSION=v0.17.2
NIMBLE_VERSION=v0.8.10

base=$(pwd)

yum install -y git curl

git clone -b $NIM_VERSION --depth 1 git://github.com/nim-lang/nim nim-$NIM_VERSION/
cd nim-$NIM_VERSION
git clone -b $NIM_VERSION --depth 1 git://github.com/nim-lang/csources csources/

cd csources
sh build.sh
cd ..
rm -rf csources
bin/nim c koch
./koch boot -d:release

export PATH=$PATH:$base/nim-$NIM_VERSION/bin/:$PATH:$base/nimble/src

cd $base

git clone -b $NIMBLE_VERSION --depth 1 git://github.com/nim-lang/nimble.git
cd nimble
nim c src/nimble
cp ./src/nimble /usr/bin/

cd $base
git clone --depth 1 git://github.com/brentp/hts-nim.git
cd hts-nim
grep -v requires hts.nimble > k.nimble && mv k.nimble hts.nimble
nimble install -y

cd $base
git clone --depth 1 git://github.com/docopt/docopt.nim.git
cd docopt.nim
nimble install -y

cd $base
git clone --depth 1 git://github.com/brentp/mosdepth.git
cd mosdepth
nim c -d:release mosdepth.nim 
cp ./mosdepth /io
