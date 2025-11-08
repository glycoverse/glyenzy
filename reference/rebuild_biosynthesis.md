# Rebuild the Biosynthetic Path of Glycans

Reconstruct the biosynthetic pathway for one or more glycans using
enzymatic reactions. This function uses a multi-target breadth-first
search to find all feasible pathways that can synthesize all the target
glycans.

## Usage

``` r
rebuild_biosynthesis(glycans, enzymes = NULL, max_steps = 20, filter = NULL)
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Can also be a single glycan. If multiple glycans are provided, the
  starting structure will be decided by the first glycan. Therefore,
  please make sure `glycans` are not of mixed glycan types.

- enzymes:

  A character vector of gene symbols, or a list of
  [`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)
  objects. If `NULL` (default), all available enzymes will be used.

- max_steps:

  Integer, maximum number of enzymatic steps to search. Default is 20.

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
gene symbols and `step` attribute indicating the step number. For
multiple targets, the graph includes all synthesis paths needed to reach
every target glycan.

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

# Rebuild the biosynthetic pathway of a single glycan
glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
path <- rebuild_biosynthesis(glycan, max_steps = 20)

# Rebuild pathways for multiple glycans
glycans <- c(
  "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
  "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
)
path <- rebuild_biosynthesis(glycans, max_steps = 20)

# View the path
igraph::as_data_frame(path, what = "edges")
#>                                                    from
#> 1                                            GalNAc(a1-
#> 2                                   Gal(b1-3)GalNAc(a1-
#> 3                       GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 4                       GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 5                       GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 6                       GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 7                       GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 8                       GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 9              Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 10             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 11             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 12             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 13             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 14             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 15             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 16             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 17  Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 18  Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 19  Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 20 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 21 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 22 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 23 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 24 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 25 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#>                                                                 to  enzyme step
#> 1                                              Gal(b1-3)GalNAc(a1- C1GALT1    1
#> 2                                  GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-  B3GNT3    2
#> 3                         Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- B4GALT1    3
#> 4                         Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- B4GALT2    3
#> 5                         Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- B4GALT3    3
#> 6                         Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- B4GALT4    3
#> 7                         Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- B4GALT5    3
#> 8                         Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- B4GALT6    3
#> 9              Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT3    4
#> 10             Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT4    4
#> 11             Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT5    4
#> 12             Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT6    4
#> 13             Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT9    4
#> 14            Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- ST3GAL3    4
#> 15            Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- ST3GAL4    4
#> 16            Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- ST3GAL6    4
#> 17 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- ST3GAL3    5
#> 18 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- ST3GAL4    5
#> 19 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1- ST3GAL6    5
#> 20 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT3    5
#> 21 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT4    5
#> 22 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT5    5
#> 23 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT6    5
#> 24 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT7    5
#> 25 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-    FUT9    5
```
