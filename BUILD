#!/usr/bin/env bash

set -e

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Everything below that is not specific to this repository lives in the family's
# tools/, reached through the BiOCamLib submodule.  What stays here is the list
# of what this repository builds, the native dependencies it vendors and the
# profiles and suites it offers -- and nothing else, because everything else was
# a copy of the same logic in every repository of the family, free to drift in
# any of them.
TOOLS="$ROOT/BiOCamLib/tools"

if [[ "${1:-}" == "README.pdf" ]]; then
  # No stylesheet of this repository's own is passed, and none exists any more:
  # the house look lives in tools/markdown.css, and a figure is sized from the
  # PNG header carried in its own data URI by tools/figbox.awk, which is what
  # README.css used to be needed for.
  bash "$TOOLS/markdown-pdf" --root "$ROOT" --title KPop
  exit 0
fi

# Release packaging and the macOS CI live in tools/release, which takes the
# project name and reads what ships from releases/MANIFEST:
#   ./BUILD package [<ver>]   assemble releases/KPop-<ver>-<os>-<arch>.tar.xz
#   ./BUILD mac-begin         tag v<CURRENT> and push it, triggering the CI
#   ./BUILD mac-end           wait for it, download the macOS binaries, package
if [[ "${1:-}" == "package" ]]; then
  bash "$TOOLS/release" package "${2:-}" --root "$ROOT" --name KPop
  exit 0
fi

if [[ "${1:-}" == "mac-begin" ]]; then
  bash "$TOOLS/release" mac-begin --root "$ROOT"
  exit 0
fi

if [[ "${1:-}" == "mac-end" ]]; then
  bash "$TOOLS/release" mac-end --root "$ROOT" --name KPop
  exit 0
fi

# Accepted targets:
#   README.pdf            Regenerate README.pdf from README.md (pandoc +
#                         headless Chrome).  Handled above; needs no build.
#   package / mac-begin / mac-end
#                         Release packaging.  Handled above; see the block
#                         just before this one.
#   <profile>             dev | dev-static | release | release-static (default dev)
#                         Builds the KPop binaries.
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
    echo "Usage: $0 [dev|dev-static|release|release-static|test|test-core|test-phylo|README.pdf|package|mac-begin|mac-end]" >&2
    exit 2
    ;;
esac

# BLAS target: arch-aware default (ARMV8 on arm64, HASWELL otherwise).  Only the
# macOS branch below consumes it; the Linux branch hardwires DYNAMIC_ARCH/PRESCOTT.
if [[ -z "${BLAS_TARGET:-}" ]]; then
  case "$(uname -m)" in
    arm64|aarch64) BLAS_TARGET="ARMV8" ;;
    *)             BLAS_TARGET="HASWELL" ;;
  esac
fi

# Always erase dune _build directory to ensure peace of mind
rm -rf _build

# ...but we want to keep our build so as not to have to recompile OpenBLAS or faiss every time
mkdir -p .build

rm -f lib/libopenblas.a
rm -f lib/libfaiss.a
rm -f lib/libinterfaiss.a

if [[ "$(uname -s)" == "Darwin" ]]; then
  # ────────────────────────────────────────────────────────────────────
  # macOS: the Linux in-binary ISA-dispatch machinery below (GNU-binutils
  # objcopy/ld -r trickery over ELF, x86 AVX multiplexing) has no Mach-O
  # analogue, and arm64 has no AVX at all.  So on macOS we build a SINGLE
  # faiss variant per architecture (generic/NEON on arm64, avx2 on Intel)
  # and a plain interfaiss shim — no run-time dispatch.  Compilers come from
  # the environment (the CI sets CC/CXX/FC to clang/clang++/gfortran); OpenMP
  # is Homebrew's libomp.  This path is exercised only on the macOS CI runner.
  # ────────────────────────────────────────────────────────────────────
  CC="${CC:-cc}"
  CXX="${CXX:-c++}"
  FC="${FC:-gfortran}"
  if command -v nproc >/dev/null 2>&1; then NPROC="$(nproc)"; else NPROC="$(sysctl -n hw.ncpu)"; fi

  # OpenBLAS: a fixed TARGET (no DYNAMIC_ARCH) for this arch, with the reference
  # LAPACK (netlib) that faiss links against.  -fno-lto sidesteps a gfortran/clang
  # LTO object mismatch on macOS; OpenMP is off inside the BLAS itself.
  if [[ -f OpenBLAS/libopenblas.a ]]; then
    cp OpenBLAS/libopenblas.a lib/
  else
    ( cd OpenBLAS
      make -j "$NPROC" libs netlib \
        CC="$CC" FC="$FC" HOSTCC="$CC" TARGET="$BLAS_TARGET" \
        NO_AVX="${NO_AVX:-0}" USE_OPENMP="${USE_OPENMP:-0}" FFLAGS="-fno-lto"
      cp libopenblas.a ../lib/ )
  fi

  # faiss: one ISA variant only.
  case "$(uname -m)" in
    arm64|aarch64) FAISS_OPT_LEVEL="generic"; FAISS_TARGET="faiss" ;;
    *)             FAISS_OPT_LEVEL="avx2";    FAISS_TARGET="faiss_avx2" ;;
  esac
  if [[ -f ".build/faiss/faiss/lib${FAISS_TARGET}.a" ]]; then
    cp ".build/faiss/faiss/lib${FAISS_TARGET}.a" lib/libfaiss.a
  else
    LIBOMP_PREFIX="$(brew --prefix libomp)"
    ( cd faiss
      cmake -D CMAKE_VERBOSE_MAKEFILE=true -D CMAKE_CXX_COMPILER="$CXX" \
        -D BLAS_LIBRARIES="$ROOT/lib/libopenblas.a" \
        -D LAPACK_LIBRARIES="$ROOT/lib/libopenblas.a" \
        -D FAISS_ENABLE_GPU=false -D FAISS_ENABLE_PYTHON=false \
        -D BUILD_TESTING=false -D CMAKE_BUILD_TYPE=Release \
        -D FAISS_OPT_LEVEL="$FAISS_OPT_LEVEL" \
        -D OpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include" \
        -D OpenMP_CXX_LIB_NAMES="omp" \
        -D OpenMP_omp_LIBRARY="${LIBOMP_PREFIX}/lib/libomp.a" \
        -B ../.build/faiss . )
    ( cd .build/faiss && make -j "$NPROC" "$FAISS_TARGET" )
    cp ".build/faiss/faiss/lib${FAISS_TARGET}.a" lib/libfaiss.a
  fi

  # interfaiss: compile the shim once (no -DIFSUF, so the entry points keep their
  # bare names — exactly what interfaiss_ocaml.c calls) and fold it together with
  # the single faiss variant into one static archive, so lib/dune's
  # `foreign_archives interfaiss` pulls faiss in too (the Linux bundle, minus the
  # ISA multiplexing).  libtool -static is the Mach-O way to merge static libs.
  ( cd lib
    LIBOMP_PREFIX="$(brew --prefix libomp)"
    "$CXX" -std=c++17 -I ../faiss/ -O3 -fPIC \
      -Xpreprocessor -fopenmp -I"${LIBOMP_PREFIX}/include" \
      -c interfaiss.cpp -o interfaiss.o
    libtool -static -o libinterfaiss.a interfaiss.o libfaiss.a
    rm -f interfaiss.o libfaiss.a )
else
  # ────────────────────────────────────────────────────────────────────
  # Linux (and any ELF/GNU-binutils host): the full machinery — OpenBLAS
  # DYNAMIC_ARCH run-time CPU dispatch, the three faiss ISA variants, and the
  # interfaiss in-binary ISA dispatch.  Left exactly as it was.
  # ────────────────────────────────────────────────────────────────────

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
  ( if [[ -f .build/faiss/faiss/libfaiss.a && -f .build/faiss/faiss/libfaiss_avx2.a && -f .build/faiss/faiss/libfaiss_avx512.a ]]; then
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
      cmake -D CMAKE_VERBOSE_MAKEFILE=true -D CMAKE_CXX_COMPILER="$(realpath ../compilers/cxx)" -D BLAS_LIBRARIES="$(realpath ../OpenBLAS/libopenblas.a)" -D FAISS_ENABLE_GPU=false -D FAISS_ENABLE_PYTHON=false -D BUILD_TESTING=false -D CMAKE_BUILD_TYPE=Release -D FAISS_OPT_LEVEL=avx512 -B ../.build/faiss .
      cd ../.build/faiss
      make -j "$(( ( $(nproc) + 1 ) / 2 ))" faiss faiss_avx2 faiss_avx512
    fi )
  cp .build/faiss/faiss/libfaiss.a        lib/libfaiss_generic.a
  cp .build/faiss/faiss/libfaiss_avx2.a   lib/libfaiss_avx2.a
  cp .build/faiss/faiss/libfaiss_avx512.a lib/libfaiss_avx512.a

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
fi

# Build everything else

# TWO VERSIONS ARE STAMPED, from two trees, because they are two different claims.
# The vendored library's is generated from the SUBMODULE's history, that being the
# tree it is compiled from, and is what `BiOCamLib.Info.info` reports when a binary
# lists what it was built against.  This repository's is generated from its own, and
# is what every binary here takes its version and date from.
#
# Each binary reports the release named in releases/CURRENT followed by this tree's
# commit-file count, so the part people cite leads and the suffix still tells two
# builds of one release apart.  Generating it is what keeps a binary from claiming a
# version its tree does not have, which is the whole reason not to write one down --
# and the six of them had drifted to six different numbers, 55, 29, 33, 48, 1 and 1,
# each with a date of its own, while shipping in one archive from one tree.
# KPop-hash2kmer is absent from the list because a name with a hyphen in it is not an
# OCaml identifier; it takes the version and date from here and keeps its own name.
bash "$TOOLS/stamp-version" --root "$ROOT/BiOCamLib" --out "$ROOT/BiOCamLib/lib/Info.ml" \
  BiOCamLib AnnoTools Cophenetic FASTools NJ Octopus Parallel RC TREx Yggdrasill
bash "$TOOLS/stamp-version" --root "$ROOT" --out "$ROOT/lib/Info.ml" --open \
  KPop KPopCount KPopCountDB KPopTwist KPopTwistDB KPopPhylo

#FLAGS="--verbose"

dune build --profile="$PROFILE" bin/KPopCount.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopCountDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwist.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwistDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopPhylo.exe $FLAGS
dune build --profile="$PROFILE" bin/KPop_hash2kmer.exe $FLAGS
dune build --profile="$PROFILE" bin/KPop_autotuner.exe $FLAGS

# When a test target is in effect we also need Yggdrasill (for the
# integration scripts) and the relevant test executables.  Those stay
# under _build/; we move only the KPop binaries to .build/.
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

mv _build/default/bin/KPopCount.exe .build/KPopCount
mv _build/default/bin/KPopCountDB.exe .build/KPopCountDB
mv _build/default/bin/KPopTwist.exe .build/KPopTwist
mv _build/default/bin/KPopTwistDB.exe .build/KPopTwistDB
mv _build/default/bin/KPopPhylo.exe .build/KPopPhylo
mv _build/default/bin/KPop_hash2kmer.exe .build/KPop-hash2kmer
mv _build/default/bin/KPop_autotuner.exe .build/KPop-autotuner

chmod 755 .build/*

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
  strip .build/KPop{Count,CountDB,Twist,TwistDB,Phylo,-hash2kmer,-autotuner}
  rm -rf _build
fi
