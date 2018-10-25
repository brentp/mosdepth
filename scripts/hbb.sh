set -e
mkdir -p built/
cp scripts/build.sh built/
docker run -t -i --rm -v `pwd`/built:/io \
    centos:centos6 \
    bash -c "yum install -y git curl wget && wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo && yum install -y devtoolset-2-gcc devtoolset-2-binutils devtoolset-2-gcc-c++ && chmod +x /io/build.sh && scl enable devtoolset-2 /io/build.sh"
echo "binary in built/"
