# Find a Biosynthesis Path Between Glycan Structures

Find biosynthetic paths from one glycan structure to another. The
default method uses known enzyme rules in a forward breadth-first
search. The virtual-enzyme method trims `to` backward to `from` and
returns every possible residue-addition order.

## Usage

``` r
path_biosynthesis(
  from,
  to,
  enzymes = NULL,
  max_steps = 10,
  filter = NULL,
  method = c("enzymatic", "virtual")
)
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
  objects. If `NULL` (default), all available enzymes will be used. Must
  be `NULL` when `method = "virtual"`.

- max_steps:

  Integer, maximum number of enzymatic steps to search. Default is 10.

- filter:

  Optional function to filter generated glycans at each step. Should
  take a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector as input and return a logical vector of the same length. For
  `method = "enzymatic"`, it filters generated products. For
  `method = "virtual"`, it filters generated precursors during backward
  trimming.

- method:

  Biosynthesis inference method. `"enzymatic"` (default) uses known
  enzyme rules in a forward search. `"virtual"` uses virtual enzymes in
  a backward search that removes terminal residues from `to` until
  `from` is reached.

## Value

An
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
object representing the synthesis path(s). Vertices represent glycan
structures with `name` attribute containing IUPAC-condensed strings.
Edges represent enzymatic reactions with `enzyme` attribute containing
gene symbols or virtual-enzyme names and `step` attribute indicating the
forward synthesis step.

## Important notes

Here are some important notes for all functions in the `glyenzy`
package.

### Applicability

Known-enzyme algorithms and enzyme information in glyenzy are applicable
only to humans, and specifically to N-glycans and O-GalNAc glycans.
Results may be inaccurate for other types of glycans (e.g., GAGs,
glycolipids) or for glycans in other species (e.g., plants, insects).

### Inclusiveness

The algorithm takes an intentionally inclusive approach, assuming that
all possible isoenzymes capable of catalyzing a given reaction may be
involved. Therefore, results should be interpreted with caution.

For example, in humans, detection of the motif "Neu5Ac(a2-3)Gal(b1-"
will return both "ST3GAL3" and "ST3GAL4". In reality, only one of them
might be active, depending on factors such as tissue specificity.

### Concrete glycans by default

Most functions only work for glycans containing **concrete** residues
(e.g., `"Glc"`, `"GalNAc"`), and not for glycans with **generic**
residues (e.g., `"Hex"`, `"HexNAc"`). Reduced-level inputs with generic
residues are supported where explicitly documented, such as
`apply_enzyme(structure_level = "basic")`,
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md),
and `path_biosynthesis()`.

### Substituents

Substituents (e.g. sulfation, phosphorylation) are not supported yet,
and the algorithms might fail for glycans with substituents. If your
glycans contain substituents, use
[`glyrepr::remove_substituents()`](https://glycoverse.github.io/glyrepr/reference/remove_substituents.html)
to get clean glycans.

### Incomplete glycan structures

If the glycan structure is incomplete or partially degraded, the result
may be misleading. Glycans with a
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
other than `"intact"` are matched with the lenient motif matching mode
in glymotif, and a warning is raised because enzyme predictions may be
less reliable.

### Starting points

For known-enzyme path inference:

- For N-glycans, the starting structure is assumed to be
  "Glc(3)Man(9)GlcNAc(2)", the N-glycan precursor transferred to Asn by
  OST.

- For O-GalNAc glycans, the starting structure is assumed to be
  "GalNAc(a1-".

- For O-GlcNAc glycans, the starting structure is assumed to be
  "GlcNAc(b1-".

- For O-Man glycans, the starting structure is assumed to be "Man(a1-".

- For O-Fuc glycans, the starting structure is assumed to be "Fuc(a1-".

- For O-Glc glycans, the starting structure is assumed to be "Glc(b1-".

## Virtual enzymes

With `method = "virtual"`, each edge is named for the residue added by
that step. Intact glycans include the linkage anomer and acceptor
position, so a beta-1,4-linked GlcNAc is labeled `"b4GlcNAcT"`. Partial
and topological glycans omit linkage information and use `"GlcNAcT"`;
basic glycans use the generic residue name, such as `"HexNAcT"`.

Virtual tracing starts N-glycans at the N-glycan core and all other
glycans at their reducing-end root residue. In `path_biosynthesis()`,
the explicit `from` glycan is always the virtual starting structure.
These networks do not apply organism-specific substrate rules and
represent structural possibilities rather than biological feasibility.

Basic structures do not retain glycan-class metadata. A basic structure
matching the generic N-glycan-core topology is therefore assumed to be
an N-glycan; use `path_biosynthesis()` with an explicit `from` when that
topology belongs to another glycan class.

## Examples

``` r
library(glyrepr)
library(glyparse)

# Find shortest path
from <- "Gal(b1-4)GlcNAc(b1-"
to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
path <- path_biosynthesis(from, to, enzymes = "ST6GAL1", max_steps = 3)

# Ignore known enzyme specificity and infer additions by trimming backward
virtual_path <- path_biosynthesis(
  "Gal(b1-3)GalNAc(a1-",
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
  method = "virtual"
)

# View the path
igraph::as_data_frame(path, what = "edges")
#>                  from                              to  enzyme step
#> 1 Gal(b1-4)GlcNAc(b1- Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1- ST6GAL1    1
```
