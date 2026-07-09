# Match Residues Added by an Enzyme

This function finds residues in glycans that match the product motifs of
a glycosyltransferase and returns their node indices.

## Usage

``` r
match_enzyme(glycans, enzyme, method = c("motif", "path"))
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector.

- enzyme:

  A glycosyltransferase
  [`enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/enzyme.md)
  or a gene symbol for one. Glycoside hydrolases are not supported.

- method:

  Method used to decide whether the enzyme is involved. `"motif"`
  matches product motifs directly in each glycan. `"path"` matches
  substrates and products from
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  results back to each glycan, which is more accurate but slower.

## Value

A list of integer vectors with the same length as `glycans`. Each
integer vector contains node indices for residues added by `enzyme` in
the corresponding glycan.

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

### Concrete glycans by default

Most functions only work for glycans containing **concrete** residues
(e.g., `"Glc"`, `"GalNAc"`), and not for glycans with **generic**
residues (e.g., `"Hex"`, `"HexNAc"`). Reduced-level inputs with generic
residues are supported where explicitly documented, such as
`apply_enzyme(structure_level = "basic")`,
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md),
and
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md).

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
glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
match_enzyme(glycan, "ST3GAL3")
#> [[1]]
#> [1] 1
#> 
match_enzyme(glycan, "ST3GAL3", method = "path")
#> [[1]]
#> [1] 1
#> 
```
