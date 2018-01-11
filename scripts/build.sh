set -euo pipefail

base=$(pwd)

yum install -y git

NIM_VERSION=v0.17.2

git clone -b $NIM_VERSION --depth 1 git://github.com/nim-lang/nim nim-$NIM_VERSION/
cd nim-$NIM_VERSION
git clone --depth 1 git://github.com/nim-lang/csources csources/

cd csources
sh build.sh
cd ..
rm -rf csources
bin/nim c koch
./koch boot -d:release

export PATH=$PATH:$base/nim-$NIM_VERSION/bin/:$PATH:$base/nimble/src

cd $base
git clone --depth 1 https://github.com/nim-lang/nimble.git
cd nimble
nim c src/nimble
src/nimble install -y

cd $base
git clone --depth 1 https://github.com/brentp/hts-nim
cd hts-nim
grep -v requires hts.nimble > k.nimble && mv k.nimble hts.nimble
nimble install -y

git clone --depth 1 https://github.com/brentp/mosdepth
cd mosdepth
nimble install -y docopt

nim c -d:release mosdepth.nim 
cp ./mosdepth /io
