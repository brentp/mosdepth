set -euo pipefail

NIM_VERSION=v0.17.2

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
#/koch nimble

export PATH=$PATH:$base/nim-$NIM_VERSION/bin/:$PATH:$base/nimble/src

cd $base

#git clone --depth 1 https://github.com/nim-lang/nimble.git
#cd nimble
curl -Lo master.tar.gz  https://github.com/nim-lang/nimble/archive/master.tar.gz
tar xzvf master.tar.gz
rm -f master.tar.gz
cd nimble-master
nim c src/nimble
cp ./src/nimble /usr/bin/
#src/nimble install -y

cd $base
curl -Lo master.tar.gz  https://github.com/brentp/hts-nim/archive/master.tar.gz
#git clone --depth 1 https://github.com/brentp/hts-nim
tar xzvf master.tar.gz
rm -rf master.tar.gz
cd hts-nim-master
grep -v requires hts.nimble > k.nimble && mv k.nimble hts.nimble
nimble install -y

cd $base
curl -Lo master.tar.gz  https://github.com/docopt/docopt.nim/archive/master.tar.gz
tar xzvf master.tar.gz
rm -rf master.tar.gz
cd docopt.nim-master
nimble install -y

#git clone --depth 1 https://github.com/brentp/mosdepth
cd $base
curl -Lo master.tar.gz  https://github.com/brentp/mosdepth/archive/master.tar.gz
tar xzvf master.tar.gz
rm -rf master.tar.gz
cd mosdepth-master
nim c -d:release mosdepth.nim 
cp ./mosdepth /io
