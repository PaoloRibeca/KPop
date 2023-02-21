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

#FLAGS="--verbose"

dune build --profile="$PROFILE" bin/KPopCount.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopCountDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwist_.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwistDB.exe $FLAGS

mv _build/default/bin/KPopCount.exe build/KPopCount
mv _build/default/bin/KPopCountDB.exe build/KPopCountDB
mv _build/default/bin/KPopTwist_.exe build/KPopTwist_
mv _build/default/bin/KPopTwistDB.exe build/KPopTwistDB

chmod 755 build/*

if [[ "$PROFILE" == "release" ]]; then
  strip build/*
  rm -rf _build
fi

cp src/KPop* build

