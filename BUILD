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

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ──────────────────────────────────────────────────────────────────────
# Release packaging targets (mirroring BioTattoo):
#   ./BUILD package [<ver>]   assemble releases/KPop-<ver>-<os>-<arch>.tar.xz
#                             from releases/MANIFEST.  <ver> is N.N.N; omitted
#                             it is read from releases/CURRENT, given it replaces
#                             CURRENT (old one kept as CURRENT~).  Builds the
#                             binaries first if none are present (release-static
#                             on Linux, release on macOS).  Set PACKAGE_PLATFORM
#                             + PACKAGE_BINDIR to package downloaded binaries for
#                             another platform (see mac-end).
#   ./BUILD mac-begin         tag v<CURRENT> on HEAD (if absent) + push it,
#                             triggering the macOS CI (build-binaries.yml).
#   ./BUILD mac-end           wait for that CI run, download the macOS binaries
#                             and 'package' them into releases/ (needs gh).
# ──────────────────────────────────────────────────────────────────────

# Assemble a distributable package (tarball) under releases/.  The package dir
# is KPop-<version>-<os>-<arch> (os/arch from uname, so the same target works on
# Linux and macOS), populated from releases/MANIFEST — a two-column list mapping
# each source path in the tree (first column) to its name inside the package
# (second column) — then tar'd and xz-compressed.  A '.build/<name>' entry is a
# compiled binary (taken from the local build tree, or from PACKAGE_BINDIR when
# packaging downloaded binaries for another platform); every other entry is
# copied from the source tree.
build_package() {
  local version_arg="${1:-}"
  local manifest="$ROOT/releases/MANIFEST"
  local current="$ROOT/releases/CURRENT"
  local version_re='^[0-9]+\.[0-9]+\.[0-9]+$'
  local version

  command -v xz >/dev/null 2>&1 \
    || { echo "(BUILD package): xz not found" >&2; exit 1; }
  [[ -f "$manifest" ]] \
    || { echo "(BUILD package): releases/MANIFEST not found" >&2; exit 1; }

  if [[ -n "$version_arg" ]]; then
    # Version given on the command line: validate it, then make it the new
    # CURRENT, preserving the previous one as CURRENT~.
    [[ "$version_arg" =~ $version_re ]] \
      || { echo "(BUILD package): '$version_arg' is not a valid version (expected N.N.N)" >&2; exit 1; }
    version="$version_arg"
    if [[ -f "$current" ]]; then mv -f "$current" "$current~"; fi
    printf '%s\n' "$version" > "$current"
  else
    # No version given: take it from CURRENT, which must exist and be valid.
    [[ -f "$current" ]] \
      || { echo "(BUILD package): no version given and releases/CURRENT not found" >&2; exit 1; }
    read -r version < "$current" || true
    [[ "$version" =~ $version_re ]] \
      || { echo "(BUILD package): releases/CURRENT does not hold a valid version (expected N.N.N)" >&2; exit 1; }
  fi

  # Target platform tag: normally the host (uname -s '-' uname -m), but a
  # cross-built package is assembled from downloaded binaries by pointing
  # PACKAGE_PLATFORM at the <os>-<arch> tag and PACKAGE_BINDIR at their
  # directory.  Dash-separated fields (e.g. KPop-1.99.1-Linux-x86_64) let the
  # parts be recovered by splitting on '-': neither os nor arch carries a dash on
  # the platforms we target, and arch keeps its underscore (x86_64).
  local platform name dir
  if [[ -n "${PACKAGE_PLATFORM:-}" ]]; then
    [[ "$PACKAGE_PLATFORM" =~ ^[A-Za-z0-9_-]+$ ]] \
      || { echo "(BUILD package): PACKAGE_PLATFORM '$PACKAGE_PLATFORM' has unexpected characters (want e.g. Darwin-arm64)" >&2; exit 1; }
    platform="$PACKAGE_PLATFORM"
  else
    platform="$(uname -s)-$(uname -m)"
  fi
  name="KPop-${version}-${platform}"
  dir="$ROOT/releases/$name"

  # Where the compiled binaries (the .build/<name> MANIFEST entries) come from.
  local bindir
  if [[ -n "${PACKAGE_BINDIR:-}" ]]; then
    # Cross-built package: take the binaries from the downloaded directory and
    # compile nothing locally.
    bindir="$PACKAGE_BINDIR"
    [[ -d "$bindir" ]] \
      || { echo "(BUILD package): PACKAGE_BINDIR '$bindir' is not a directory" >&2; exit 1; }
  else
    # Native package: build the binaries first if none are present.  release-static
    # bundles libgmp/musl so the Linux binaries run on hosts lacking them; macOS
    # has no static libSystem, so -ccopt -static won't link there — use plain
    # release.
    bindir="$ROOT/.build"
    local profile
    case "$(uname -s)" in
      Darwin) profile="release" ;;
      *)      profile="release-static" ;;
    esac
    if [[ ! -x "$bindir/KPopCount" ]]; then
      echo "(BUILD package): .build/KPopCount missing — building with '$profile' first ..."
      ( cd "$ROOT" && bash BUILD "$profile" )
    fi
  fi

  echo "(BUILD package): assembling $name ..."
  rm -rf "$dir"
  mkdir -p "$dir"

  # Populate from MANIFEST: first column = source path in the tree, second
  # column = name inside the package.  A missing source fails the build loudly
  # rather than yielding a half-empty archive.
  local src dst from
  while read -r src dst; do
    [[ -z "$src" ]] && continue          # blank line
    case "$src" in \#*) continue;; esac  # comment line
    [[ -n "$dst" ]] \
      || { echo "(BUILD package): malformed MANIFEST line for '$src' (need two columns)" >&2; exit 1; }
    case "$src" in
      .build/*) from="$bindir/${src#.build/}" ;;
      *)       from="$ROOT/$src" ;;
    esac
    [[ -f "$from" ]] \
      || { echo "(BUILD package): '$src' -> '$from' not found — build it first?" >&2; exit 1; }
    cp -p "$from" "$dir/$dst"
    # Restore the executable bit on binaries: GitHub artifacts are zipped and
    # lose Unix permissions, so a downloaded binary may arrive non-executable.
    case "$src" in .build/*) chmod +x "$dir/$dst" ;; esac
  done < "$manifest"

  # tar the directory, then compress the archive with xz.  Run from within
  # releases/ so the archive holds paths relative to it (no releases/ prefix).
  echo "(BUILD package): creating $name.tar.xz ..."
  ( cd "$ROOT/releases" && tar cf "$name.tar" "$name" && xz -f "$name.tar" )
  echo "(BUILD package): done — releases/$name.tar.xz"
}

# Kick off the macOS CI build.  Reads the version from releases/CURRENT, forms
# the tag v<version>, creates it on HEAD if it does not exist yet, and pushes it
# to origin — which triggers .github/workflows/build-binaries.yml.  Pair with
# 'bash BUILD mac-end' to collect and package the resulting binaries.
build_mac_begin() {
  local current="$ROOT/releases/CURRENT"
  local version_re='^[0-9]+\.[0-9]+\.[0-9]+$'
  local version tag
  [[ -f "$current" ]] \
    || { echo "(BUILD mac-begin): releases/CURRENT not found" >&2; exit 1; }
  read -r version < "$current" || true
  [[ "$version" =~ $version_re ]] \
    || { echo "(BUILD mac-begin): releases/CURRENT does not hold a valid version (expected N.N.N)" >&2; exit 1; }
  tag="v$version"
  if git -C "$ROOT" rev-parse -q --verify "refs/tags/$tag" >/dev/null; then
    echo "(BUILD mac-begin): tag $tag already exists — not recreating."
  else
    echo "(BUILD mac-begin): creating tag $tag on HEAD ..."
    git -C "$ROOT" tag "$tag"
  fi
  echo "(BUILD mac-begin): pushing $tag to origin (triggers the macOS CI) ..."
  git -C "$ROOT" push origin "$tag"
  echo "(BUILD mac-begin): done — run 'bash BUILD mac-end' to wait for the build,"
  echo "                   then download and package the macOS binaries."
}

# Collect the macOS CI build and package it.  Waits for the build-binaries run
# of tag v<version> (version from releases/CURRENT) to finish, downloads the
# Darwin-arm64 and Darwin-x86_64 binaries, and runs 'package' for each — so the
# macOS archives land in releases/ next to the Linux one.  Needs an authenticated
# gh (GitHub CLI): run 'gh auth login' first.
build_mac_end() {
  local current="$ROOT/releases/CURRENT"
  local version_re='^[0-9]+\.[0-9]+\.[0-9]+$'
  local version tag
  command -v gh >/dev/null 2>&1 \
    || { echo "(BUILD mac-end): gh (GitHub CLI) not found — needed to fetch the CI artifacts" >&2; exit 1; }
  gh auth status >/dev/null 2>&1 \
    || { echo "(BUILD mac-end): gh is not authenticated — run 'gh auth login' first" >&2; exit 1; }
  [[ -f "$current" ]] \
    || { echo "(BUILD mac-end): releases/CURRENT not found" >&2; exit 1; }
  read -r version < "$current" || true
  [[ "$version" =~ $version_re ]] \
    || { echo "(BUILD mac-end): releases/CURRENT does not hold a valid version (expected N.N.N)" >&2; exit 1; }
  tag="v$version"

  # Find the CI run for this tag.  It may take a few seconds to register after
  # the tag push, so poll briefly before giving up.
  echo "(BUILD mac-end): locating the macOS CI run for $tag ..."
  local run_id="" i
  for ((i = 0; i < 30; i++)); do
    run_id="$(gh run list --workflow build-binaries.yml --event push --branch "$tag" \
                --limit 1 --json databaseId --jq '.[0].databaseId // empty' 2>/dev/null || true)"
    [[ -n "$run_id" ]] && break
    sleep 5
  done
  [[ -n "$run_id" ]] \
    || { echo "(BUILD mac-end): no CI run found for $tag — did 'bash BUILD mac-begin' run?" >&2; exit 1; }

  # Wait for it to finish; --exit-status makes gh fail if the run concluded badly.
  echo "(BUILD mac-end): waiting for run $run_id to finish ..."
  gh run watch "$run_id" --exit-status

  # Download each macOS artifact into its own scratch directory and package it.
  # 'package' with no version argument reads releases/CURRENT, so it never
  # re-shuffles CURRENT and the archives share the Linux build's version.
  local staging
  staging="$(mktemp -d)"
  local arch platform
  for arch in arm64 x86_64; do
    platform="Darwin-$arch"
    echo "(BUILD mac-end): downloading and packaging $platform ..."
    mkdir -p "$staging/$platform"
    gh run download "$run_id" --name "KPop-$platform" --dir "$staging/$platform"
    PACKAGE_PLATFORM="$platform" PACKAGE_BINDIR="$staging/$platform" \
      bash "$ROOT/BUILD" package
  done
  rm -rf "$staging"
  echo "(BUILD mac-end): done — macOS packages are in releases/."
}

if [[ "${1:-}" == "package" ]]; then
  build_package "${2:-}"
  exit 0
fi

if [[ "${1:-}" == "mac-begin" ]]; then
  build_mac_begin
  exit 0
fi

if [[ "${1:-}" == "mac-end" ]]; then
  build_mac_end
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

# Emit version info for both BiOCamLib and KPop.  The date is formatted by git
# itself (--date=format) rather than `date -d @<ts>`, which is GNU-only and fails
# on macOS's BSD date; the version stays the git file-change count.
cd BiOCamLib && echo -e "include (\n  struct\n    let info = {\n      Tools.Argv.name = \"BiOCamLib\";\n      version = \"$(git log --pretty=format: --name-only | awk '{if ($0!="") print}' | wc -l)\";\n      date = \"$(git log -1 --format=%ad --date=format:'%d-%b-%Y')\"\n    }\n  end\n)" > lib/Info.ml && cd ..
echo -e "include (\n  struct\n    let info = {\n      BiOCamLib.Tools.Argv.name = \"KPop\";\n      version = \"$(git log --pretty=format: --name-only | awk '{if ($0!="") print}' | wc -l)\";\n      date = \"$(git log -1 --format=%ad --date=format:'%d-%b-%Y')\"\n    }\n  end\n)" > lib/Info.ml

#FLAGS="--verbose"

dune build --profile="$PROFILE" bin/KPopCount.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopCountDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwist.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopTwistDB.exe $FLAGS
dune build --profile="$PROFILE" bin/KPopPhylo.exe $FLAGS
dune build --profile="$PROFILE" bin/KPop_hash2kmer.exe $FLAGS

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
  strip .build/KPop{Count,CountDB,Twist,TwistDB,Phylo,-hash2kmer}
  rm -rf _build
fi
