#!/bin/sh

# You might need to change this first line depending on your installation
rm -f BiOCamLib
ln -s ../BiOCamLib

dune build bin/KPopCount.exe
dune build bin/KPopCountDB.exe
dune build bin/kPopTwist.exe
dune build bin/KPopTwistDB.exe

chmod 755 _build/default/bin/KPopCount.exe _build/default/bin/KPopCountDB.exe _build/default/bin/kPopTwist.exe _build/default/bin/KPopTwistDB.exe

rm -f KPopCount KPopCountDB kPopTwist KPopTwistDB
ln -s _build/default/bin/KPopCount.exe KPopCount
ln -s _build/default/bin/KPopCountDB.exe KPopCountDB
ln -s _build/default/bin/kPopTwist.exe kPopTwist
ln -s _build/default/bin/KPopTwistDB.exe KPopTwistDB

