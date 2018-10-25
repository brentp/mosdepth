set -euo pipefail

NIM_VERSION=v0.19.0

base=$(pwd)

git clone -b $NIM_VERSION --depth 1 git://github.com/nim-lang/nim nim-$NIM_VERSION/
cd nim-$NIM_VERSION

sh build_all.sh

export PATH=$PATH:$base/nim-$NIM_VERSION/bin/

cd $base

cd $base
git clone --depth 1 git://github.com/brentp/mosdepth.git
cd mosdepth
nimble install -y
nim c -d:release mosdepth.nim 
cp ./mosdepth /io
