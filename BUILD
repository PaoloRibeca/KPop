#!/bin/sh

#####################################################################
# You might need to change this line depending on your installation #
BIOCAMLIB_PATH=../BiOCamLib
#####################################################################

rm -f BiOCamLib
ln -s "$BIOCAMLIB_PATH" BiOCamLib

PROFILE="$1"
if [ -z "$PROFILE" ]; then
  PROFILE="dev"
fi

dune build --profile="$PROFILE" bin/KPopCount.exe
dune build --profile="$PROFILE" bin/KPopCountDB.exe
dune build --profile="$PROFILE" bin/KPopTwist_.exe
dune build --profile="$PROFILE" bin/KPopTwistDB.exe

chmod 755 _build/default/bin/KPopCount.exe _build/default/bin/KPopCountDB.exe _build/default/bin/KPopTwist_.exe _build/default/bin/KPopTwistDB.exe
if [ "$PROFILE" = "release" ]; then
  strip _build/default/bin/KPopCount.exe _build/default/bin/KPopCountDB.exe _build/default/bin/KPopTwist_.exe _build/default/bin/KPopTwistDB.exe
fi

rm -f KPopCount KPopCountDB KPopTwist_ KPopTwistDB
ln -s _build/default/bin/KPopCount.exe KPopCount
ln -s _build/default/bin/KPopCountDB.exe KPopCountDB
ln -s _build/default/bin/KPopTwist_.exe KPopTwist_
ln -s _build/default/bin/KPopTwistDB.exe KPopTwistDB

