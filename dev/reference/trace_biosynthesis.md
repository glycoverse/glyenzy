# Trace the Biosynthetic Path of Glycans

Reconstruct biosynthetic pathways for one or more glycans. The default
method uses known enzyme rules in a forward, multi-target breadth-first
search. To infer structure-driven paths without enzyme specificity, use
[`trace_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis_virtual.md).

## Usage

``` r
trace_biosynthesis(
  glycans,
  enzymes = NULL,
  max_steps = NULL,
  filter = NULL,
  max_virtual_steps = 0L
)
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
  objects. If `NULL` (default), all available enzymes will be used.

- max_steps:

  Integer, maximum number of enzymatic steps to search, or `NULL` to
  infer it from the largest target glycan. The inferred value is the
  maximum monosaccharide and substituent count plus 4 for N-glycans, and
  the maximum count minus 1 for other glycans.

- filter:

  Optional function to filter generated glycans at each step. Should
  take a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector as input and return a logical vector of the same length. It
  filters generated products.

- max_virtual_steps:

  Integer, maximum number of target-directed virtual enzyme steps
  allowed when no fully enzymatic path exists. Default is `0L`, which
  disables virtual fallback. See the "Virtual fallback" section for more
  details.

## Value

A `glyenzy_biosynthesis_network` object inheriting from
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
and representing the synthesis path(s). Vertices represent glycan
structures, with IUPAC-condensed strings in the `name` attribute and a
logical `target` attribute indicating whether each vertex is a target
glycan. At most one directed edge connects each substrate and product.
Every edge has an integer `step`, a logical `is_virtual`, a scalar
display label in `enzyme`, and a list-valued `enzymes` attribute
containing concrete enzyme candidates. For known reactions, `enzyme`
combines the candidates with `" / "`; for virtual reactions, it contains
the structural virtual-enzyme name and `enzymes` is empty unless
concrete candidates were annotated.

For multiple targets, the graph includes all synthesis paths needed to
reach every target glycan.

## Input compatibility

`trace_biosynthesis()`,
[`trace_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis_virtual.md),
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md),
and
[`path_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis_virtual.md)
accept only glycans that share one monosaccharide type (`"concrete"` or
`"generic"`) and one structure level (`"intact"` or `"topological"`).
Structures with mixed generic and concrete residues, partial structures,
and missing structures are not supported. In the `path_*()` functions,
`from` and `to` must have the same monosaccharide type and structure
level.

Use
[`glyrepr::convert_to_generic()`](https://glycoverse.github.io/glyrepr/reference/convert_to_generic.html)
to standardize monosaccharide types or
[`glyrepr::remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.html)
to standardize structures at the topological level.

## Virtual fallback

Sometimes the biosynthesis network of a glycan cannot be fully resolved;
i.e., some enzymatic steps are not inferred to be catalyzed by any known
enzyme ("bad" steps). By default, an error is raised for these glycans.

`max_virtual_steps` provides a fallback for these glycans. For a "bad"
step, a virtual enzyme is assigned to allow the algorithm to continue.
For example, for the O-GalNAc core 5 "GalNAc(a1-3)GalNAc(a1-", an
"a3GalNAcT" is assigned to the step that adds the a3 GalNAc. Unsupported
sulfate additions similarly use `"3SulfoT"`, `"6SulfoT"`, or
`"?SulfoT"`.

Therefore, `max_virtual_steps` can also be interpreted as "the maximum
number of glycosidic bonds or sulfate transfers that cannot be assigned
by a known enzyme." Increasing this number loosens the criteria.

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
residues (e.g., `"Hex"`, `"HexNAc"`). Inputs with generic residues are
supported where explicitly documented, such as `trace_biosynthesis()`
and
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md).

### Substituents

Sulfate substituents are supported. Other substituents, such as
phosphorylation and methylation, are not supported. Use
[`glyrepr::remove_substituents()`](https://glycoverse.github.io/glyrepr/reference/remove_substituents.html)
when unsupported substituents are present.

### Incomplete or non-concrete glycan structures

Some functions accept glycan structures that are incomplete, partially
degraded, or contain generic residues. When accepted, glycans with a
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

- For GlcCer glycans, the starting structure is assumed to be "Glc(b1-",

- For GalCer glycans, the starting structure is assumed to be "Gal(b1-"

## Examples

``` r
library(glyrepr)
library(glyparse)

# Rebuild the biosynthetic pathway of a single glycan
glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
path <- trace_biosynthesis(glycan)

# Rebuild pathways for multiple glycans
glycans <- c(
  "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
  "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
)
path <- trace_biosynthesis(glycans)

# View the path
igraph::as_data_frame(path, what = "edges")
#>                                                   from
#> 1                                           GalNAc(a1-
#> 2                                  Gal(b1-3)GalNAc(a1-
#> 3                      GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 4             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 5             Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 6 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#>                                                                to
#> 1                                             Gal(b1-3)GalNAc(a1-
#> 2                                 GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 3                        Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 4             Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 5            Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> 6 Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#>                                            enzyme is_virtual step      enzymes
#> 1                                         C1GALT1      FALSE    1      C1GALT1
#> 2                                          B3GNT3      FALSE    2       B3GNT3
#> 3 B4GALT1 / B4GALT2 / B4GALT3 / B4GALT4 / B4GALT5      FALSE    3 B4GALT1,....
#> 4                FUT3 / FUT4 / FUT5 / FUT6 / FUT9      FALSE    4 FUT3, FU....
#> 5                               ST3GAL4 / ST3GAL6      FALSE    4 ST3GAL4,....
#> 6                              FUT5 / FUT6 / FUT7      FALSE    5 FUT5, FU....
```
