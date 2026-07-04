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
  # Render to a temp PDF; only move it into place once it is a valid PDF,
  # so a Chrome crash never leaves a stale README.pdf looking like success.
  PDF="$(mktemp --suffix=.pdf)"
  LOG="$(mktemp)"
  trap 'rm -rf "$HTML" "$UDD" "$PDF" "$LOG"' EXIT
  # --embed-resources (pandoc >= 2.19) supersedes the older --self-contained;
  # use whichever this pandoc advertises so the build also works on older installs.
  EMBED=--embed-resources
  pandoc --help 2>/dev/null | grep -q -- --embed-resources || EMBED=--self-contained
  pandoc README.md -f gfm -t html5 --standalone "$EMBED" --mathml \
         --css README.css --metadata title="KPop" -o "$HTML"
  # --disable-background-networking is essential: the GCM/SSL phone-home it
  # otherwise attempts stalls — and on this Chrome crashes (dangling raw_ptr)
  # — the headless render of a large image-heavy page.  The rest are headless
  # hygiene so the render is fully self-contained.
  if "$CHROME" --headless=new --no-sandbox --disable-gpu \
       --disable-dev-shm-usage --disable-background-networking \
       --disable-default-apps --disable-extensions --disable-sync \
       --disable-component-update --no-first-run --metrics-recording-only \
       --user-data-dir="$UDD" --no-pdf-header-footer \
       --print-to-pdf="$PDF" "$HTML" 2>"$LOG" \
     && [[ -s "$PDF" ]] && [[ "$(head -c4 "$PDF")" == "%PDF" ]]; then
    mv "$PDF" README.pdf
    echo "BUILD: wrote README.pdf"
  else
    echo "BUILD: README.pdf FAILED — pandoc or headless Chrome error:" >&2
    sed 's/^/    /' "$LOG" >&2
    exit 1
  fi
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

# Build OpenBLAS.  DYNAMIC_ARCH bundles every x86-64 kernel and selects one at
# run time from CPUID (Nehalem..SapphireRapids/Zen), while a PRESCOTT baseline
# keeps the common code portable -- so the binaries run on old/VM CPU models
# without AVX yet still reach AVX-512 speed on modern hardware.
( if [[ -f OpenBLAS/libopenblas.a ]]; then
    cp OpenBLAS/libopenblas.a lib/
  else
    cd OpenBLAS
    make -j "$(nproc)" CC="$(realpath ../compilers/cc)" FC="$(realpath ../compilers/fortran)" HOSTCC="$(realpath ../compilers/cc)" DYNAMIC_ARCH=1 DYNAMIC_OLDER=1 TARGET=PRESCOTT CROSS=1 NO_SHARED=1
    cp libopenblas.a ../lib/
  fi )

# Build faiss in three ISA variants (generic + avx2 + avx512).  FAISS itself does
# not dispatch inside a single static library (its Python packaging selects a
# per-ISA shared object at import), so we build all three and dispatch ourselves
# below.  OPT_LEVEL=avx512 un-excludes the generic and avx2 targets too.
( if [[ -f build/faiss/faiss/libfaiss.a && -f build/faiss/faiss/libfaiss_avx2.a && -f build/faiss/faiss/libfaiss_avx512.a ]]; then
    :
  else
    # faiss requires cmake >= 3.24; fail early and clearly if the cmake on PATH is
    # too old (providing a suitable cmake is the user's responsibility, not ours)
    v=$(cmake --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+' | head -1)
    if [[ -z "$v" ]] || (( 10#${v%%.*} * 100 + 10#${v##*.} < 324 )); then
      echo "BUILD: faiss requires cmake >= 3.24 (found '${v:-none}' on PATH)" >&2
      exit 1
    fi
    cd faiss
    cmake -D CMAKE_VERBOSE_MAKEFILE=true -D CMAKE_CXX_COMPILER="$(realpath ../compilers/cxx)" -D BLAS_LIBRARIES="$(realpath ../OpenBLAS/libopenblas.a)" -D FAISS_ENABLE_GPU=false -D FAISS_ENABLE_PYTHON=false -D BUILD_TESTING=false -D CMAKE_BUILD_TYPE=Release -D FAISS_OPT_LEVEL=avx512 -B ../build/faiss .
    cd ../build/faiss
    # avx512 translation units are memory-heavy; we halve the number of workers
    # to cap peak RAM
    make -j "$(( ( $(nproc) + 1 ) / 2 ))" faiss faiss_avx2 faiss_avx512
  fi )
cp build/faiss/faiss/libfaiss.a        lib/libfaiss_generic.a
cp build/faiss/faiss/libfaiss_avx2.a   lib/libfaiss_avx2.a
cp build/faiss/faiss/libfaiss_avx512.a lib/libfaiss_avx512.a

# Build interfaiss.  Compile the FAISS shim once per ISA (entry points suffixed
# via -include interfaiss_variant.h) and bundle each with its matching FAISS variant
# into one relocatable object, so the three coexist and are picked at run time by
# interfaiss_dispatch.c.  Two things make this safe in a single static binary:
#   1. each bundle's ISA CODE (strong .text functions) is localised, keeping only the
#      7 suffixed entry points global, so the variants don't clash;
#   2. the specialised (avx2/avx512) variants would otherwise SIGILL at startup, because
#      C++ static initialisers run unconditionally before main and carry their variant's
#      ISA.  So we DATA-share: the variants' strong data globals are weakened (the generic
#      copy wins) and their .init_array/.fini_array are stripped -- only the SSE2 generic
#      constructors run, initialising the one shared set of globals that the dispatched
#      kernels then use.  The generic path stays AVX2-free, so it runs on any x86-64.
( cd lib
  CXX="$(realpath ../compilers/cxx)"
  CC="$(realpath ../compilers/cc)"
  ENTRIES="interfaiss_create_flat_index interfaiss_create_PQ_index interfaiss_create_HNSW_index interfaiss_query_index interfaiss_add_data_to_index interfaiss_train_index interfaiss_free_index"
  for v in generic avx2 avx512; do
    case "$v" in
      generic) ISA="" ;;
      avx2)    ISA="-mavx2 -mfma -mf16c -mpopcnt" ;;
      avx512)  ISA="-mavx2 -mfma -mf16c -mavx512f -mavx512cd -mavx512vl -mavx512dq -mavx512bw -mpopcnt" ;;
    esac
    "$CXX" -I ../faiss/ -O3 -fPIC -fopenmp $ISA -include interfaiss_variant.h -DIFSUF=_$v -c interfaiss.cpp -o if_$v.o
    ld -r -o bundle_$v.o if_$v.o libfaiss_$v.a
    : > keep_$v.txt; for e in $ENTRIES; do echo "${e}_$v" >> keep_$v.txt; done
    nm --defined-only -g bundle_$v.o | awk '$2 == "T" {print $NF}' | grep -vxFf keep_$v.txt > localize_$v.txt
    objcopy --localize-symbols=localize_$v.txt bundle_$v.o bundle_${v}_loc.o
    if [ "$v" != generic ]; then
      nm --defined-only -g bundle_${v}_loc.o | awk '$2 ~ /^[BDGRS]$/ {print $NF}' > weaken_$v.txt
      objcopy --weaken-symbols=weaken_$v.txt bundle_${v}_loc.o
      objcopy --wildcard --remove-section '.init_array*' --remove-section '.fini_array*' --remove-section '.ctors*' --remove-section '.dtors*' bundle_${v}_loc.o
    fi
  done
  "$CC" -I ../faiss/ -O3 -fPIC -c interfaiss_dispatch.c -o interfaiss_dispatch.o
  rm -f libinterfaiss.a
  ar rcs libinterfaiss.a interfaiss_dispatch.o bundle_generic_loc.o bundle_avx2_loc.o bundle_avx512_loc.o
  rm -f if_*.o bundle_*.o keep_*.txt localize_*.txt weaken_*.txt interfaiss_dispatch.o libfaiss_generic.a libfaiss_avx2.a libfaiss_avx512.a )

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
dune build --profile="$PROFILE" bin/KPop_hash2kmer.exe $FLAGS

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
mv _build/default/bin/KPop_hash2kmer.exe build/KPop-hash2kmer

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
  strip build/KPop{Count,CountDB,Twist,TwistDB,Phylo,-hash2kmer}
  rm -rf _build
fi

