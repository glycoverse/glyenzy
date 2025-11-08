# Find Synthesis Path Between Glycan Structures

Find a synthesis path from one glycan structure to another using
enzymatic reactions. This function uses breadth-first search to find the
shortest path or all possible paths within a given number of steps.

## Usage

``` r
find_synthesis_path(from, to, enzymes = NULL, max_steps = 10, filter = NULL)
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
  [`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)
  objects. If `NULL` (default), all available enzymes will be used.

- max_steps:

  Integer, maximum number of enzymatic steps to search. Default is 10.

- filter:

  Optional function to filter generated glycans at each step. Should
  take a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector as input and return a logical vector of the same length. It
  will be applied to all the generated glycans at each BFS step for
  pruning.

## Value

An
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
object representing the synthesis path(s). Vertices represent glycan
structures with `name` attribute containing IUPAC-condensed strings.
Edges represent enzymatic reactions with `enzyme` attribute containing
gene symbols and `step` attribute indicating the step number.

## Important notes

Here are some important notes for all functions in the `glyenzy`
package.

### Applicability

All algorithms and enzyme information in glyenzy are applicable only to
humans, and specifically to N-glycans and O-GalNAc glycans. Results may
be inaccurate for other types of glycans (e.g., GAGs, glycolipids) or
for glycans in other species (e.g., plants, insects).

### Inclusiveness

The algorithm takes an intentionally inclusive approach, assuming that
all possible isoenzymes capable of catalyzing a given reaction may be
involved. Therefore, results should be interpreted with caution.

For example, in humans, detection of the motif "Neu5Ac(a2-3)Gal(b1-"
will return both "ST3GAL3" and "ST3GAL4". In reality, only one of them
might be active, depending on factors such as tissue specificity.

### Only "concrete" glycans

The function only works for glycans containing **concrete** residues
(e.g., `"Glc"`, `"GalNAc"`), and not for glycans with **generic**
residues (e.g., `"Hex"`, `"HexNAc"`).

### Substituents

Subtituents (e.g. sulfation, phosphorylation) is not supported yet, and
the algorithms might fail for glycans with subtituents. If your glycans
contains substituents, use
[`glyrepr::remove_substituents()`](https://glycoverse.github.io/glyrepr/reference/remove_substituents.html)
to get clean glycans.

### Incomplete glycan structures

If the glycan structure is incomplete or partially degraded, the result
may be misleading.

### Starting points

- For N-glycans, the starting structure is assumed to be
  "Glc(3)Man(9)GlcNAc(2)", the N-glycan precursor transfered to Asn by
  OST.

- For O-GalNAc glycans, the starting structure is assumed to be
  "GalNAc(a1-".

- For O-GlcNAc glycans, the starting structure is assumed to be
  "GlcNAc(b1-".

- For O-Man glycans, the starting structure is assumed to be "Man(a1-".

- For O-Fuc glycans, the starting structure is assumed to be "Fuc(a1-".

- For O-Glc glycans, the starting structure is assumed to be "Glc(b1-".

## Examples

``` r
library(glyrepr)
library(glyparse)

# Find shortest path
from <- "Gal(b1-4)GlcNAc(b1-"
to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
path <- find_synthesis_path(from, to, enzymes = "ST6GAL1", max_steps = 3)

# View the path
igraph::as_data_frame(path, what = "edges")
#>                  from                              to  enzyme step
#> 1 Gal(b1-4)GlcNAc(b1- Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1- ST6GAL1    1
```
