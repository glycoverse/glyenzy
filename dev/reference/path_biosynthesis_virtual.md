# Find a Virtual Biosynthesis Path Between Glycan Structures

Infer every structure-driven biosynthetic path from `from` to `to` by
trimming `to` backward to `from`. Unlike
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md),
this does not require known enzyme rules.

## Usage

``` r
path_biosynthesis_virtual(from, to, enzymes = NULL, annotate_enzymes = FALSE)
```

## Arguments

- from:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  scalar, or a character string supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  The starting glycan structure.

- to:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  scalar, or a character string supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  The target glycan structure.

- enzymes:

  A character vector of gene symbols, or a list of
  [`enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/enzyme.md)
  objects. Used only when `annotate_enzymes` is `TRUE`; if `NULL`, all
  available enzymes are considered.

- annotate_enzymes:

  Whether to annotate each virtual transition with concrete enzymes
  whose rules can perform it. Defaults to `FALSE`.

## Value

An
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
object representing the synthesis path(s). Vertices contain
IUPAC-condensed strings in `name`; edges have a forward `step` and
virtual-enzyme `enzyme` attribute. When `annotate_enzymes` is `TRUE`,
`concrete_enzymes` is a list of character vectors containing every
candidate concrete enzyme for each transition.

## Virtual enzymes

Each edge is named for the residue added by that step. Intact glycans
include the linkage anomer and acceptor position, so a beta-1,4-linked
GlcNAc is labeled `"b4GlcNAcT"`. Partial and topological glycans omit
linkage information and use `"GlcNAcT"`; basic glycans use the generic
residue name, such as `"HexNAcT"`.

Virtual tracing starts N-glycans at the N-glycan core and all other
glycans at their reducing-end root residue. In
`path_biosynthesis_virtual()`, the explicit `from` glycan is always the
virtual starting structure. These networks represent structural
possibilities rather than biological feasibility.

Basic structures do not retain glycan-class metadata. A basic structure
matching the generic N-glycan-core topology is therefore assumed to be
an N-glycan; use `path_biosynthesis_virtual()` with an explicit `from`
when that topology belongs to another glycan class.

## Examples

``` r
virtual_path <- path_biosynthesis_virtual(
  "Gal(b1-3)GalNAc(a1-",
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
)

annotated_path <- path_biosynthesis_virtual(
  "Gal(b1-3)GalNAc(a1-",
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
  annotate_enzymes = TRUE
)
```
