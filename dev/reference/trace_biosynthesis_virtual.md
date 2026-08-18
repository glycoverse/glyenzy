# Trace a Virtual Biosynthetic Path of Glycans

Reconstruct every structure-driven biosynthetic path for one or more
glycans by trimming terminal residues and sulfate groups backward.
Unlike
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md),
this does not require known enzyme rules.

## Usage

``` r
trace_biosynthesis_virtual(glycans, enzymes = NULL, annotate_enzymes = FALSE)
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
  [`enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/enzyme.md)
  objects. Used only when `annotate_enzymes` is `TRUE`; if `NULL`, all
  available enzymes are considered.

- annotate_enzymes:

  Whether to annotate each virtual transition with concrete enzymes
  whose rules can perform it. Defaults to `FALSE`.

## Value

A `glyenzy_virtual_biosynthesis_network` object inheriting from
`glyenzy_biosynthesis_network` and
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html).
Vertices contain IUPAC-condensed strings in `name` and a logical
`target` attribute indicating whether each vertex is a target glycan;
edges have a forward `step` and virtual-enzyme `enzyme` attribute. When
`annotate_enzymes` is `TRUE`, `concrete_enzymes` is a list of character
vectors containing every candidate concrete enzyme for each transition.

## Virtual enzymes

Each edge is named for the residue added by that step. Intact glycans
include the linkage anomer and acceptor position, so a beta-1,4-linked
GlcNAc is labeled `"b4GlcNAcT"`. Partial and topological glycans omit
linkage information and use `"GlcNAcT"`. Generic or mixed topological
glycans use their preserved generic residue names, such as `"HexNAcT"`.

Sulfation is represented as its own atomic transition. Sulfate additions
at positions 3 and 6 use `"3SulfoT"` and `"6SulfoT"`; an unknown or
other position uses `"?SulfoT"`. A sulfated terminal residue is
therefore desulfated before the residue itself can be trimmed.

Virtual tracing starts N-glycans at the N-glycan core and all other
glycans at their reducing-end root residue. Sulfates are removed from
these automatically selected starts. In
[`path_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis_virtual.md),
the explicit `from` glycan is always the virtual starting structure,
including any sulfate groups it contains; those sulfates must also occur
in `to`. These networks represent structural possibilities rather than
biological feasibility.

Generic structures do not retain glycan-class metadata. A generic
structure matching the N-glycan-core topology is therefore assumed to be
an N-glycan; use
[`path_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis_virtual.md)
with an explicit `from` when that topology belongs to another glycan
class.

## Important notes

Here are some important notes for all functions in the `glyenzy`
package.

### Applicability

Known-enzyme algorithms and enzyme information in glyenzy are applicable
only to humans. Curated coverage is strongest for N-glycans and
O-glycans and also includes selected glycosphingolipid headgroups and
other glycan contexts. Lipid and protein aglycones are not represented,
so glycolipid rules model the carbohydrate headgroup with ceramide
omitted. Results may be inaccurate for unsupported glycan contexts or
other species (e.g., plants, insects).

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
residues (e.g., `"Hex"`, `"HexNAc"`). Inputs with generic or mixed
residues are supported where explicitly documented, such as
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
and
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md).

### Substituents

Sulfate substituents are supported. Other substituents, such as
phosphorylation and methylation, are not supported. Use
[`glyrepr::remove_substituents()`](https://glycoverse.github.io/glyrepr/reference/remove_substituents.html)
when unsupported substituents are present.

### Incomplete or non-concrete glycan structures

If the glycan structure is incomplete, partially degraded, or contains
generic or mixed residues, the result may be misleading. Glycans with a
[`glyrepr::get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.html)
other than `"intact"`, or with a
[`glyrepr::get_mono_type()`](https://glycoverse.github.io/glyrepr/reference/get_mono_type.html)
other than `"concrete"`, are matched with the lenient motif matching
mode in glymotif. A warning is raised because enzyme predictions may be
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

virtual_path <- trace_biosynthesis_virtual(
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
)

annotated_path <- trace_biosynthesis_virtual(
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
  annotate_enzymes = TRUE
)
```
