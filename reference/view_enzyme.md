# View Residues Added by an Enzyme

Visualize where an enzyme contributes residues to a glycan structure.

## Usage

``` r
view_enzyme(glycan, enzyme)
```

## Arguments

- glycan:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  or a glycan structure string supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- enzyme:

  A glycosyltransferase
  [`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)
  or a gene symbol for one. Glycoside hydrolases are not supported.

## Value

A `ggplot` object returned by
[`glydraw::draw_cartoon()`](https://glycoverse.github.io/glydraw/reference/draw_cartoon.html).
If no match is found, the glycan is drawn without highlighted residues
and a cli alert is emitted.

## Details

`view_enzyme()` matches one `enzyme` against one `glycan` with the same
matching rules used by
[`match_enzyme()`](https://glycoverse.github.io/glyenzy/reference/match_enzyme.md),
then draws the glycan with the matched residues highlighted.

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

Substituents (e.g. sulfation, phosphorylation) are not supported yet,
and the algorithms might fail for glycans with substituents. If your
glycans contain substituents, use
[`glyrepr::remove_substituents()`](https://glycoverse.github.io/glyrepr/reference/remove_substituents.html)
to get clean glycans.

### Incomplete glycan structures

If the glycan structure is incomplete or partially degraded, the result
may be misleading.

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

## See also

[`match_enzyme()`](https://glycoverse.github.io/glyenzy/reference/match_enzyme.md),
[`glydraw::draw_cartoon()`](https://glycoverse.github.io/glydraw/reference/draw_cartoon.html)

## Examples

``` r
glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")

if (FALSE) { # \dontrun{
view_enzyme(glycan, "ST3GAL3")
} # }
```
