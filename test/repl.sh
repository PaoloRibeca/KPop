#!/usr/bin/env bash
#
# repl.sh [<script.ml> ...]
#
# An OCaml toplevel with this repository's libraries already linked in, for
# asking `lib/` a question without writing a test to ask it with.  Given no
# argument it is interactive; given files it runs them and exits.
#
#     bash test/repl.sh                     # a prompt, with KPop open
#     bash test/repl.sh /tmp/probe.ml       # run a script and stop
#     echo 'Space.Distance.of_string "euclidean";;' | bash test/repl.sh
#
# WHY IT IS A CUSTOM TOPLEVEL AND NOT `ocaml`.  The switches this project builds
# in are musl-linked and static, so the stock bytecode toplevel cannot dlopen
# zarith's stubs and `#require` fails on the first C dependency.  `ocamlmktop
# -custom` links the stubs and the C libraries into the toplevel itself, which
# is the same thing `bin/dune` does for the binaries -- and here that is a
# longer list than elsewhere in the family: OpenBLAS, faiss and the interfaiss
# C++ shim, with gfortran, quadmath, stdc++ and OpenMP under them.
#
# AND WHY `-ccopt -static` ON TOP OF `-custom`.  The two settle different halves:
# `-custom` archives the OCaml stubs into the executable and leaves the C
# libraries under them to be resolved at load time, which on this musl runtime
# libstdc++ and libgomp are not.  Static brings them in too.
#
# WHAT IT IS FOR: a question about a library function whose answer you want
# before deciding what the test should assert -- what `Clustering.valley_radius`
# makes of a set of embeddings, what a metric does to a distance.  It is NOT a
# substitute for the suites: nothing here is checked by anything, so an answer
# found this way belongs in a test before it is relied on.
#
# The bytecode archives are asked for under the RELEASE profile, because `dev`
# turns warnings into errors and lib/SparseNJ.ml carries two declarations it
# does not use -- which stops the toplevel building for a reason that has
# nothing to do with the toplevel.
#
# The toplevel is cached beside the build and rebuilt when either library is
# newer than it, so the usual run costs nothing.

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
B="$ROOT/_build/default"
TOP="$ROOT/.build/KPop-repl"
BIO_OBJS="$B/BiOCamLib/lib/.BiOCamLib.objs/byte"
KPOP_OBJS="$B/lib/.KPop.objs/byte"
BIO_CMA="$B/BiOCamLib/lib/BiOCamLib.cma"
KPOP_CMA="$B/lib/KPop.cma"

# The libraries are built native for the binaries, so the bytecode archives the
# toplevel needs are asked for by name rather than assumed to be there.
( cd "$ROOT" && dune build --profile release lib/KPop.cma BiOCamLib/lib/BiOCamLib.cma ) || {
  echo "repl.sh: could not build the bytecode libraries" >&2; exit 1
}

if [[ ! -x "$TOP" || "$KPOP_CMA" -nt "$TOP" || "$BIO_CMA" -nt "$TOP" ]]; then
  mkdir -p "$(dirname "$TOP")"
  # Both directories are searched: _build/default/lib holds the stubs archive
  # dune compiles from lib/*.c, and lib/ holds the two vendored ones BUILD puts
  # there.  There is no -lfaiss: faiss is archived INSIDE libinterfaiss.a rather
  # than beside it, which is why the binaries name only the two.  They are given
  # in dependency order -- interfaiss over OpenBLAS -- because a static archive
  # is searched only for symbols already wanted when the linker reaches it.
  ocamlfind ocamlmktop -custom -ccopt -static \
    -package str,unix,zarith -linkpkg \
    -I "$BIO_OBJS" -I "$KPOP_OBJS" "$BIO_CMA" "$KPOP_CMA" \
    -cclib -L"$B/lib" -cclib -L"$ROOT/lib" \
    -cclib -linterfaiss -cclib -lopenblas \
    -cclib -lstdc++ -cclib -lgfortran -cclib -lquadmath -cclib -lgomp \
    -o "$TOP" 2> >(grep -v 'warning: Using' >&2) || {
      echo "repl.sh: could not build the toplevel" >&2; exit 1
    }
fi

# `open KPop` up front, so that a script says `Clustering.valley_radius` and not
# `KPop.Clustering.valley_radius`, which is how the library reads everywhere else.
PRELUDE="$(mktemp)"
trap 'rm -f "$PRELUDE"' EXIT
# Terminated, because this same file is fed to the toplevel on stdin when a
# script is given and an unterminated phrase there swallows what follows it.
cat > "$PRELUDE" <<'EOF'
open BiOCamLib;;
open Better;;
open KPop;;
EOF

if [[ $# -eq 0 ]]; then
  exec "$TOP" -I "$BIO_OBJS" -I "$KPOP_OBJS" -init "$PRELUDE"
fi

# A SCRIPT IS `#use`d AND NOT CONCATENATED.  The toplevel ignores `-init` once it
# is given a file to run, so the opens have to reach the script some other way --
# and pasting them in front of it would shift every line number in the file,
# which is the one thing a person debugging a script needs to be able to trust.
# `#use` reports errors against the file's own lines.
for f in "$@"; do
  [[ -r "$f" ]] || { echo "repl.sh: cannot read '$f'" >&2; exit 1; }
  printf '#use "%s";;\n' "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" >> "$PRELUDE"
done
exec "$TOP" -I "$BIO_OBJS" -I "$KPOP_OBJS" -noprompt < "$PRELUDE"
