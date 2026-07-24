# Identify Potentially Involved Enzymes

This function returns all possible isoenzymes associated with the
biosynthetic steps of the input glycan. Note that this function ignores
the residues in glycans that cannot be matched to any enzyme rules.

## Usage

``` r
find_enzyme(glycans, return_list = NULL, method = c("motif", "path"))
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- return_list:

  If `NULL` (default), return a list of character vectors when `glycans`
  has length greater than 1, and a single character vector when
  `glycans` has length 1. Set to `TRUE` to always return a list. This
  can be useful when you are working programmatically with unknown input
  length. Note that when `return_list = FALSE` and
  `length(glycans) > 1`, an error will be thrown.

- method:

  Method used to infer enzyme involvement. `"motif"` checks product
  motifs directly in each glycan. `"path"` extracts enzymes from
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/reference/trace_biosynthesis.md)
  results, which is more accurate but slower.

## Value

A character vector or a list of character vectors (see `return_list`
parameter), each containing the names of enzymes involved in the
biosynthesis of the corresponding glycan.

## Important notes

Here are some important notes for all functions in the `glyenzy`
package.

### Applicability

Known-enzyme algorithms and enzyme information in glyenzy are applicable
only to humans, and specifically to N-glycans and O-glycans. Results may
be inaccurate for other types of glycans (e.g., GAGs, glycolipids) or
for glycans in other species (e.g., plants, insects).

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
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/reference/trace_biosynthesis.md),
and
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/reference/path_biosynthesis.md).

### Substituents

Sulfate substituents are supported. Other substituents, such as
phosphorylation and methylation, are not supported. Use
[`glyrepr::remove_substituents()`](https://glycoverse.github.io/glyrepr/reference/remove_substituents.html)
when unsupported substituents are present.

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

## Examples

``` r
library(glyrepr)
library(glyparse)

# Use `glycan_structure()`
glycans <- auto_parse(c(
  "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
))
find_enzyme(glycans)
#> [[1]]
#>  [1] "ALG13"  "ALG14"  "ALG1"   "ALG2"   "ALG11"  "ALG3"   "ALG9"   "ALG12" 
#>  [9] "ALG6"   "ALG8"   "ALG10"  "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2"
#> [17] "MAN1C1" "MAN2A1" "MAN2A2" "MGAT1"  "DPAGT1"
#> 
#> [[2]]
#>  [1] "ALG13"  "ALG14"  "ALG1"   "ALG2"   "ALG11"  "ALG3"   "ALG9"   "ALG12" 
#>  [9] "ALG6"   "ALG8"   "ALG10"  "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2"
#> [17] "MAN1C1" "MAN2A1" "MAN2A2" "MGAT1"  "MGAT2"  "DPAGT1"
#> 

# Or use characters directly
find_enzyme("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
#>  [1] "ALG13"  "ALG14"  "ALG1"   "ALG2"   "ALG11"  "ALG3"   "ALG9"   "ALG12" 
#>  [9] "ALG6"   "ALG8"   "ALG10"  "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2"
#> [17] "MAN1C1" "MAN2A1" "MAN2A2" "MGAT1"  "DPAGT1"

# Use reconstructed biosynthesis paths
find_enzyme(glycans, method = "path")
#> [[1]]
#>  [1] "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2" "MAN1C1" "MGAT1"  "MAN2A1"
#>  [9] "MAN2A2" "ALG13"  "ALG14"  "ALG1"   "ALG2"   "ALG11"  "ALG3"   "ALG9"  
#> [17] "ALG12"  "ALG6"   "ALG8"   "ALG10" 
#> 
#> [[2]]
#>  [1] "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2" "MAN1C1" "MGAT1"  "MAN2A1"
#>  [9] "MAN2A2" "MGAT2"  "ALG13"  "ALG14"  "ALG1"   "ALG2"   "ALG11"  "ALG3"  
#> [17] "ALG9"   "ALG12"  "ALG6"   "ALG8"   "ALG10" 
#> 
```
