dune build test/PELT.exe
dune build bin/Parallel.exe
dune build bin/KPopCount.exe
dune build bin/KPopCountDB.exe
dune build bin/kPopTwist.exe
dune build bin/KPopTwistDB.exe
dune build bin/CovidDB.exe

chmod u+w _build/default/bin/{Parallel.exe,KPopCount.exe,KPopCountDB.exe,kPopTwist.exe,KPopTwistDB.exe,CovidDB.exe} _build/default/test/PELT.exe

rm PELT Parallel KPopCount KPopCountDB kPopTwist KPopTwistDB CovidDB
ln -s _build/default/test/PELT.exe PELT
ln -s _build/default/bin/Parallel.exe Parallel
ln -s _build/default/bin/KPopCount.exe KPopCount
ln -s _build/default/bin/KPopCountDB.exe KPopCountDB
ln -s _build/default/bin/kPopTwist.exe kPopTwist
ln -s _build/default/bin/KPopTwistDB.exe KPopTwistDB
ln -s _build/default/bin/CovidDB.exe CovidDB

