set -e
mkdir -p built/
cp scripts/build.sh built/
docker run -t -i --rm -v `pwd`/built:/io \
    phusion/holy-build-box-64:latest /hbb_exe/activate-exec \
    bash /io/build.sh
echo "binary in built/"
