#!/usr/bin/env bash

set -e

PROFILE="$1"
if [[ "$PROFILE" == "" ]]; then
  PROFILE="dev"
fi

# Always erase build directory to ensure peace of mind
rm -rf _build

rm -rf build
mkdir build

# Emit version info for both BiOCamLib and KPop
cd BiOCamLib && echo -e "include (\n  struct\n    let info = {\n      Tools.Argv.name = \"BiOCamLib\";\n      version = \"$(git log --pretty=format: --name-only | awk '{if ($0!="") print}' | wc -l)\";\n      date = \"$(date -d "@$(git log -1 --format="%at")" +%d-%b-%Y)\"\n    }\n  end\n)" > lib/Info.ml && cd ..
echo -e "include (\n  struct\n    let info = {\n      BiOCamLib.Tools.Argv.name = \"KPop\";\n      version = \"$(git log --pretty=format: --name-only | awk '{if ($0!="") print}' | wc -l)\";\n      date = \"$(date -d "@$(git log -1 --format="%at")" +%d-%b-%Y)\"\n    }\n  end\n)" > lib/Info.ml

#FLAGS="--verbose"

dune build --profile="$PROFILE" bin/KPopCount.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopCountDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwist_.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwistDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopPhylo.exe $FLAGS

mv _build/default/bin/KPopCount.exe build/KPopCount
mv _build/default/bin/KPopCountDB.exe build/KPopCountDB
mv _build/default/bin/KPopTwist_.exe build/KPopTwist_
mv _build/default/bin/KPopTwistDB.exe build/KPopTwistDB
mv _build/default/bin/KPopPhylo.exe build/KPopPhylo

chmod 755 build/*

if [[ "$PROFILE" == "release" ]]; then
  strip build/*
  rm -rf _build
fi

cp src/KPop* build

