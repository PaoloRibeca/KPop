(lang dune 3.0)

(name KPop)

(generate_opam_files true)

(source
 (github PaoloRibeca/KPop))

(authors "Paolo Ribeca")

(maintainers "Paolo Ribeca")

(license LICENSE)

(documentation https://github.com/PaoloRibeca/KPop#readme)

(package
 (name KPop)
 (synopsis "Unleash the full power of your k-mers!")
 (description "A set of tools to unleash the full power of your k-mers")
 (depends ocaml dune)
 (tags
  (topics k-mers bioinformatics statistics OCaml)))

; See the complete stanza docs at https://dune.readthedocs.io/en/stable/dune-files.html#dune-project
