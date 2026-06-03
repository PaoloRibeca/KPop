#!/usr/bin/env bash

set -e

# ──────────────────────────────────────────────────────────────────────
# Special target: regenerate README.pdf from README.md
#   ./BUILD README.pdf
# Markdown -> self-contained HTML (pandoc, GitHub-flavoured, image + CSS
# embedded) -> PDF (headless Chrome/Chromium, the same engine GitHub's
# "print to PDF" uses).  Needs: pandoc and google-chrome/chromium.
# ──────────────────────────────────────────────────────────────────────
if [[ "${1:-}" == "README.pdf" ]]; then
  cd "$(dirname "${BASH_SOURCE[0]}")"
  command -v pandoc >/dev/null 2>&1 || { echo "BUILD: pandoc not found" >&2; exit 1; }
  CHROME=""
  for c in google-chrome google-chrome-stable chromium chromium-browser; do
    if command -v "$c" >/dev/null 2>&1; then CHROME="$c"; break; fi
  done
  [[ -n "$CHROME" ]] || { echo "BUILD: no google-chrome/chromium found" >&2; exit 1; }
  HTML="$(mktemp --suffix=.html)"
  # Private profile dir so this Chrome never blocks on the default-profile
  # singleton lock held by an unrelated Chrome already running on the host.
  UDD="$(mktemp -d)"
  trap 'rm -rf "$HTML" "$UDD"' EXIT
  # --embed-resources (pandoc >= 2.19) supersedes the older --self-contained;
  # use whichever this pandoc advertises so the build also works on older installs.
  EMBED=--embed-resources
  pandoc --help 2>/dev/null | grep -q -- --embed-resources || EMBED=--self-contained
  pandoc README.md -f gfm -t html5 --standalone "$EMBED" --mathml \
         --css README.css --metadata title="KPop" -o "$HTML"
  "$CHROME" --headless=new --no-sandbox --disable-gpu --no-pdf-header-footer \
            --user-data-dir="$UDD" \
            --print-to-pdf=README.pdf "$HTML" 2>/dev/null
  echo "BUILD: wrote README.pdf"
  exit 0
fi

# Accepted targets:
#   README.pdf            Regenerate README.pdf from README.md (pandoc +
#                         headless Chrome).  Handled above; needs no build.
#   <profile>             dev | dev-static | release | release-static (default dev)
#                         Builds the four KPop binaries.
#   test                  Builds binaries + Yggdrasill + every test executable,
#                         then runs them all.  Fails on any test failure.
#   test-core             Builds the four binaries, then runs
#                         test/integration_core.sh which replays the Quick
#                         Start tutorial both without and with k-mer
#                         selection.  Catches regressions in the
#                         KPopCount -> KPopCountDB -> KPopTwist -> KPopTwistDB
#                         core pipeline.
#   test-phylo            Builds binaries + test/Phylo.exe, then runs the
#                         OCaml phylo-invariants exe and the
#                         test/integration_splits.sh shell suite.
case "$1" in
  "" | "dev" | "dev-static" | "release" | "release-static")
    PROFILE="${1:-dev}"
    DO_TESTS=""
    ;;
  "test")
    PROFILE="release-static"
    DO_TESTS="all"
    ;;
  "test-core")
    PROFILE="release-static"
    DO_TESTS="core"
    ;;
  "test-phylo")
    PROFILE="release-static"
    DO_TESTS="phylo"
    ;;
  *)
    echo "Usage: $0 [dev|dev-static|release|release-static|test|test-core|test-phylo|README.pdf]" >&2
    exit 2
    ;;
esac

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
dune build --profile="$PROFILE" bin/KPopPhylo.exe $FLAGS

# When a test target is in effect we also need Yggdrasill (for the
# integration scripts) and the relevant test executables.  Those stay
# under _build/; we move only the four KPop binaries to build/.
if [[ -n "$DO_TESTS" ]]; then
  case "$DO_TESTS" in
    "all")
      dune build --profile="$PROFILE" BiOCamLib/bin/Yggdrasill.exe $FLAGS
      dune build --profile="$PROFILE" \
        test/CA.exe test/RSVD.exe test/Epsilon.exe \
        test/Cluster.exe test/Phylo.exe $FLAGS
      ;;
    "phylo")
      dune build --profile="$PROFILE" test/Phylo.exe $FLAGS
      ;;
    "core")
      # No Yggdrasill or test exes needed; the four KPop binaries suffice
      :
      ;;
  esac
fi

mv _build/default/bin/KPopCount.exe build/KPopCount
mv _build/default/bin/KPopCountDB.exe build/KPopCountDB
mv _build/default/bin/KPopTwist.exe build/KPopTwist
mv _build/default/bin/KPopTwistDB.exe build/KPopTwistDB
mv _build/default/bin/KPopPhylo.exe build/KPopPhylo

chmod 755 build/*

# Run tests before any _build cleanup.  set -e propagates failures.
if [[ -n "$DO_TESTS" ]]; then
  echo
  echo "=========================================================="
  echo "  Running tests (target '$1')"
  echo "=========================================================="
  case "$DO_TESTS" in
    "phylo")
      bash test/integration_build.sh
      _build/default/test/Phylo.exe
      echo
      bash test/integration_splits.sh
      ;;
    "core")
      bash test/integration_core.sh
      ;;
    "all")
      bash test/integration_build.sh
      _build/default/test/CA.exe
      # RSVD.exe / Epsilon.exe / Cluster.exe are diagnostic tools rather
      # than self-checking unit tests: invoking them on the fixtures and
      # checking for a clean exit is the smoke test.  Stdout (large
      # comparison / diagnostic tables) is discarded; stderr surfaces
      # any crash.  RSVD's -d 10 keeps dims + oversampling well below
      # n_samples; Epsilon's -o 1 skips the O(n^2) Part 2 table
      _build/default/test/RSVD.exe -d 10 test/Primer/Train-5 >/dev/null
      _build/default/test/Epsilon.exe -o 1 test/Primer/Classes-5 >/dev/null
      _build/default/test/Cluster.exe -1 test/Primer/Classes-5 >/dev/null
      _build/default/test/Phylo.exe
      echo
      bash test/integration_core.sh
      echo
      bash test/integration_splits.sh
      ;;
  esac
fi

# Stripping + _build cleanup is for release flows that don't need _build kept
# around for testing.  Test targets leave _build intact for repeat runs.
if [[ -z "$DO_TESTS" ]] && \
   [[ "$PROFILE" == "release" || "$PROFILE" == "release-static" ]]; then
  strip build/{KPopCount,KPopCountDB,KPopTwist,KPopTwistDB}
  rm -rf _build
fi

