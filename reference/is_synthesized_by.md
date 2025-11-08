# Determine Whether a Glycan Is Synthesized by a Given Enzyme

Glycans are produced through a series of enzymatic reactions. This
function checks whether a specific enzyme participates in the
biosynthesis of a given glycan (or glycans).

## Usage

``` r
is_synthesized_by(glycans, enzyme)
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- enzyme:

  An
  [`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)
  or a gene symbol.

## Value

A logical vector of the same length as `glycans`.

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

## Algorithm

The basic approach is straightforward: for each reaction rule associated
with the enzyme, the function checks whether the corresponding product
motif appears in the glycan. If any rule matches, the function returns
`TRUE`.

For N-glycans, additional logic is applied to handle special cases.
Products of **MGAT1** are often further trimmed by glycoside hydrolases,
meaning that the final glycan product may no longer contain the original
motif. In these cases, the function instead looks for specific motif
markers to determine enzyme involvement.

## Examples

``` r
library(glyrepr)
library(glyparse)

# Use `glycan_structure()` and `enzyme()`
glycan <- auto_parse("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
is_synthesized_by(glycan, enzyme("ST6GAL1"))
#> [1] TRUE

# Or use characters directly
is_synthesized_by("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-", "ST6GAL1")
#> [1] TRUE

# Vectorized input
glycans <- c(
  "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
  "Gal(b1-4)GlcNAc(b1-"
)
is_synthesized_by(glycans, "ST6GAL1")
#> [1]  TRUE FALSE
```
