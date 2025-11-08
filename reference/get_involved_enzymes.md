# Identify Potentially Involved Enzymes

This function returns all possible isoenzymes associated with the
biosynthetic steps of the input glycan. Note that this function ignores
the residues in glycans that cannot be matched to any enzyme rules.

## Usage

``` r
get_involved_enzymes(glycans, return_list = NULL)
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

## Value

A character vector or a list of character vectors (see `return_list`
parameter), each containing the names of enzymes involved in the
biosynthesis of the corresponding glycan.

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

# Use `glycan_structure()`
glycans <- auto_parse(c(
  "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
))
get_involved_enzymes(glycans)
#> [[1]]
#> [1] "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2" "MAN1C1" "MAN2A1" "MAN2A2"
#> [9] "MGAT1" 
#> 
#> [[2]]
#>  [1] "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2" "MAN1C1" "MAN2A1" "MAN2A2"
#>  [9] "MGAT1"  "MGAT2" 
#> 

# Or use characters directly
get_involved_enzymes("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
#> [1] "MOGS"   "GANAB"  "MAN1B1" "MAN1A1" "MAN1A2" "MAN1C1" "MAN2A1" "MAN2A2"
#> [9] "MGAT1" 
```
