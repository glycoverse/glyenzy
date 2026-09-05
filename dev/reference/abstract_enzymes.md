# Abstract enzymes for N-glycan biosynthesis

Return a separate collection of simplified N-glycan reaction rules.
These enzymes represent reaction activities rather than individual
genes. They are never included in
[`db_enzymes()`](https://glycoverse.github.io/glyenzy/dev/reference/db_enzymes.md),
[`enzymes_from_rnaseq()`](https://glycoverse.github.io/glyenzy/dev/reference/enzymes_from_rnaseq.md),
or the default enzyme sets of other functions. Individual activities can
also be retrieved by name with
[`enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/enzyme.md),
for example `enzyme("ManII")`.

## Usage

``` r
abstract_enzymes(return_str = FALSE)
```

## Arguments

- return_str:

  If `FALSE` (default), return a named list of enzyme objects. If
  `TRUE`, return their activity names, which are accepted by
  [`enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/enzyme.md)
  and by `enzyme` or `enzymes` arguments. An `enzymes` collection can
  use activity names or objects, but cannot mix abstract and concrete
  enzymes.

## Value

A named list of `glyenzy_abstract_enzyme` objects, or a character vector
when `return_str = TRUE`. Each object also inherits from
`glyenzy_gt_enzyme` or `glyenzy_gh_enzyme`, and `glyenzy_enzyme`.

## Details

The collection contains ManI, GnTI, ManII, GnTII, GnTIII, a6FucT, GnTIV,
GnTV, b4GalT, iGnT, a3SiaT, and GlcH: 12 enzymes with 32 reaction rules.
All have `glycan_type = "N"` and `species = "unspecified"`. The
additional `localization` field records `"cis"`, `"medial"`, `"trans"`,
or `"ER"`; it is descriptive and does not impose compartment order
during tracing.

The trimming activities follow the corresponding concrete enzymes in the
ordinary database. ManI combines the 16 distinct whole-glycan rules of
MAN1B1, MAN1A1, MAN1A2, and MAN1C1. ManII uses the shared MAN2A1/MAN2A2
rule, including its core alignment and rejection of bisecting GlcNAc or
a galactosylated alpha1-3-arm GlcNAc. It removes either the alpha1-3 or
alpha1-6 mannose first, preserving both Man4 intermediates. These rules
give the same trimming reactions and network topology as their concrete
counterparts; enzyme labels remain abstract.

Core transferase rules also follow the concrete database, including
their acceptors, alignment, and rejects: GnTI uses MGAT1; GnTII uses
MGAT2; GnTIII uses MGAT3; a6FucT uses FUT8; GnTIV uses MGAT4A/MGAT4B;
and GnTV uses the N-glycan rule shared by MGAT5/MGAT5B. The
O-mannose-specific MGAT5B rule is not included. GnTIII adds bisecting
GlcNAc and has `localization = "medial"` in this abstract model.

FUT8, MGAT1, and MGAT2 have no explicit rejects in the concrete
database; their acceptor and alignment restrictions are retained instead
of adding new exclusions. GnTIII rejects galactosylated alpha1-3-arm
antennae. GnTIV and GnTV retain their bisecting-GlcNAc and
antenna-specific rejects. A reject excludes the matching reaction site,
not every site on a glycan.

iGnT rejects the first LacNAc extension on both the beta1-2 and beta1-4
antennae of the alpha1-3 arm, without excluding sites on another arm.
Because only iGnT initiates LacNAc extension in this collection, that
arm remains unextended in networks generated from the default precursor.
The B3GNT4 reject for terminal Gal-GlcNAc-Gal-GlcNAc additionally
prevents a second LacNAc extension. Thus the alpha1-6 arm can extend
once, but not repeatedly through this activity.

a3SiaT retains the ST3GAL4/ST3GAL6 reject for alpha1-3 fucose on its
Gal-GlcNAc acceptor. b4GalT has no explicit reject in the corresponding
N-glycan B4GALT rules. Only restrictions for the represented N-glycan
activity are imported; unrelated O-glycan and glycolipid reactions are
excluded. These terminal activities retain the table acceptors, so this
collection is not a replacement for the complete concrete database.

Ambiguous linkages use the existing motif-matching semantics, including
potentially compatible reject matches. Product-motif statistics do not
apply substrate rejects and do not establish that a product is
reachable.

GlcH combines the three glucose-removal steps represented by MOGS and
GANAB in the ordinary database. It removes one residue per step from the
complete Glc3Man9GlcNAc2 precursor, retaining Glc2 and Glc1
intermediates. Thus the collection can trace from the usual
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
starting structure without adding ordinary enzymes or changing the
start.

Pass activities explicitly, for example `enzyme("ManI")` to
[`apply_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/apply_enzyme.md),
or activity names or the whole list from `abstract_enzymes()` to
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md).
Abstract reactions are explicit enzyme rules and have
`is_virtual = FALSE` in biosynthesis networks. Existing restrictions on
GH involvement statistics and residue matching still apply; this
collection does not add GH support to those interfaces.

## References

The simplified Golgi rules follow the user-supplied reaction table based
on the N-glycosylation model of Spahn et al. (2016), with branch
notation described by Liang et al. (2020). The supplied table takes
precedence over differences in published versions, except for the
concrete trimming and core-transferase rules, added GnTIII activity, and
terminal-activity rejects described above. See the bundled rule-data
attribution.

Spahn et al. (2016). A Markov chain model for N-linked protein
glycosylation.
[doi:10.1016/j.ymben.2015.10.007](https://doi.org/10.1016/j.ymben.2015.10.007)
.

Liang et al. (2020). A Markov model of glycosylation elucidates isozyme
specificity and glycosyltransferase interactions for glycoengineering.
[doi:10.1016/j.crbiot.2020.01.001](https://doi.org/10.1016/j.crbiot.2020.01.001)
.

## Examples

``` r
abstract_enzymes(return_str = TRUE)
#>  [1] "ManI"   "GnTI"   "ManII"  "GnTII"  "GnTIII" "a6FucT" "GnTIV"  "GnTV"  
#>  [9] "b4GalT" "iGnT"   "a3SiaT" "GlcH"  
enzyme("ManII")
#> 
#> ── Enzyme: ManII ───────────────────────────────────────────────────────────────
#> ℹ Type: "GH" (Glycoside hydrolase)
#> ℹ Species: "unspecified"
#> ℹ Glycan type: "N"
#> ℹ Abstract enzyme; localization: "medial"
#> 
#> ── Rules (1) ──
#> 
#> → Rule 1: core alignment
#> Acceptor:
#> "GlcNAc(b1-2)Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#> Product: "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#> Rejects:
#> "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#> "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

man5 <- paste0(
  "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]",
  "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
)
apply_enzyme(man5, abstract_enzymes()[["GnTI"]])
#> <glycan_structure[1]>
#> [1] GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> # Unique structures: 1
trace_biosynthesis(man5, enzymes = abstract_enzymes(), max_virtual_steps = 0L)
#> IGRAPH 48ff499 DN-- 15 19 -- 
#> + attr: name (v/c), target (v/l), enzyme (e/c), is_virtual (e/l), step
#> | (e/n), enzymes (e/x)
#> + edges from 48ff499 (vertex names):
#> [1] Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-->Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [2] Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-         ->Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-         
#> + ... omitted several edges
```
