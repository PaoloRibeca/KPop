#!/usr/bin/env bash

set -e

PROFILE="$1"
if [[ "$PROFILE" == "" ]]; then
  PROFILE="dev"
fi

if [[ "$BLAS_TARGET" == "" ]]; then
  BLAS_TARGET="HASWELL"
fi

# Always erase dune _build directory to ensure peace of mind
rm -rf _build

# ...but we want to keep our build so as not to have to recompile OpenBLAS or faiss every time
mkdir -p build

rm -f lib/libopenblas.a
rm -f lib/libfaiss.a
rm -f lib/libinterfaiss.a

# Build OpenBLAS
( if [[ -f OpenBLAS/libopenblas.a ]]; then
    cp OpenBLAS/libopenblas.a lib/
  else
    cd OpenBLAS
    make -j "$(nproc)" CC="$(realpath ../compilers/cc)" FC="$(realpath ../compilers/fortran)" HOSTCC="$(realpath ../compilers/cc)" TARGET="$BLAS_TARGET" CROSS=1
    cp libopenblas.a ../lib/
  fi )

# Build faiss
( if [[ -f build/faiss/faiss/libfaiss_avx2.a ]]; then
    cp build/faiss/faiss/libfaiss_avx2.a lib/libfaiss.a
  else
    cd faiss
    cmake -D CMAKE_VERBOSE_MAKEFILE=true -D CMAKE_CXX_COMPILER="$(realpath ../compilers/cxx)" -D BLAS_LIBRARIES="$(realpath ../OpenBLAS/libopenblas.a)" -D FAISS_ENABLE_GPU=false -D FAISS_ENABLE_PYTHON=false -D BUILD_TESTING=false -D CMAKE_BUILD_TYPE=Release -D FAISS_OPT_LEVEL=avx2 -B ../build/faiss .
    cd ../build/faiss
    make -j "$(nproc)" faiss_avx2
    cp faiss/libfaiss_avx2.a ../../lib/libfaiss.a
  fi )

# Build interfaiss
( cd lib
  ../compilers/cxx -I ../faiss/ -O3 -fPIC -fopenmp -c -o libinterfaiss.o interfaiss.cpp
  ar rcs libinterfaiss.a libinterfaiss.o
  rm -f libinterfaiss.o )

# Build everything else

# Emit version info for both BiOCamLib and KPop
cd BiOCamLib && echo -e "include (\n  struct\n    let info = {\n      Tools.Argv.name = \"BiOCamLib\";\n      version = \"$(git log --pretty=format: --name-only | awk '{if ($0!="") print}' | wc -l)\";\n      date = \"$(date -d "@$(git log -1 --format="%at")" +%d-%b-%Y)\"\n    }\n  end\n)" > lib/Info.ml && cd ..
echo -e "include (\n  struct\n    let info = {\n      BiOCamLib.Tools.Argv.name = \"KPop\";\n      version = \"$(git log --pretty=format: --name-only | awk '{if ($0!="") print}' | wc -l)\";\n      date = \"$(date -d "@$(git log -1 --format="%at")" +%d-%b-%Y)\"\n    }\n  end\n)" > lib/Info.ml

#FLAGS="--verbose"

dune build --profile="$PROFILE" bin/KPopCount.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopCountDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwist.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwistDB.exe $FLAGS

mv _build/default/bin/KPopCount.exe build/KPopCount
mv _build/default/bin/KPopCountDB.exe build/KPopCountDB
mv _build/default/bin/KPopTwist.exe build/KPopTwist
mv _build/default/bin/KPopTwistDB.exe build/KPopTwistDB

chmod 755 build/*

if [[ "$PROFILE" == "release" || "$PROFILE" == "release-static" ]]; then
  strip build/{KPopCount,KPopCountDB,KPopTwist,KPopTwistDB}
  rm -rf _build
fi

