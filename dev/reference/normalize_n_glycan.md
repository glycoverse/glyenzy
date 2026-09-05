# Normalize Topological N-Glycans

Normalize selected Gal and Neu5Ac distributions in topological
N-glycans. For the supported cases, glycans with the same numbers of
antennae, Gal, and Neu5Ac are mapped to a common representative
topology.

## Usage

``` r
normalize_n_glycan(glycans)
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

## Value

A
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector of normalized glycans.

## Why does this function exist?

[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
and related functions support topological glycans. However, N-glycan
annotations can assign residues to specific antennae even when the
available data do not resolve those positions. These assignments can
introduce unsupported distinctions into a biosynthesis network.

For example, consider the two glycans below.

Glycan 1:

    Neu5Ac - Gal - GlcNAc
                          \
    Neu5Ac - Gal - GlcNAc - Man
                               \
                                Man - GlcNAc - GlcNAc -
                               /
             Gal - GlcNAc - Man

Glycan 2:

    Neu5Ac - Gal - GlcNAc
                          \
             Gal - GlcNAc - Man
                               \
                                Man - GlcNAc - GlcNAc -
                               /
    Neu5Ac - Gal - GlcNAc - Man

Both are tri-antennary glycans with two sialic acids. When mass
spectrometry data do not resolve which antennae carry the sialic acids,
either topology may be used to annotate the same unresolved structure.
Nevertheless, tracing these annotations can produce different
biosynthetic paths and introduce redundant branches when multiple
glycans are traced together.

Glycan 3 below illustrates the problem. It has one sialic acid and could
serve as a precursor to a disialylated glycan. Of the two topologies
above, only Glycan 2 can be obtained from Glycan 3 by adding a single
sialic acid. Tracing Glycan 1 and Glycan 3 together therefore requires
separate paths rather than connecting them through a single sialylation
step. When the sialic acid positions are unresolved, this distinction
reflects the annotation rather than evidence for separate biosynthetic
routes.

Glycan 3:

             Gal - GlcNAc
                          \
             Gal - GlcNAc - Man
                               \
                                Man - GlcNAc - GlcNAc -
                               /
    Neu5Ac - Gal - GlcNAc - Man

This function addresses the problem through a normalization convention.
It applies rules to reposition Gal and Neu5Ac. In this example, mapping
Glycan 1 to Glycan 2 allows the monosialylated Glycan 3 to connect to
the disialylated representative through a single sialylation step. This
convention does not imply that distinct isomers are biologically
equivalent.

## Algorithm

The algorithm repositions only Gal and Neu5Ac residues. To explain the
algorithm, we name the branches for tri- and tetra-antennary glycans.

Tri-antennary:

    A: GlcNAc
              \
    B: GlcNAc - Man
                    \
                     Man - GlcNAc - GlcNAc -
                    /
    C: GlcNAc - Man

Tetra-antennary:

    A: GlcNAc
              \
    B: GlcNAc - Man
                    \
                     Man - GlcNAc - GlcNAc -
                    /
    C: GlcNAc - Man
              /
    D: GlcNAc

First, Gal residues are repositioned in the following cases:

- Three antennae, one Gal: positioned on C.

- Three antennae, two Gal: positioned on A and C.

- Four antennae, two Gal: positioned on A and C.

Next, Neu5Ac residues are repositioned in the following cases:

- Three antennae, two Gal, one Neu5Ac: positioned on the C Gal.

- Three antennae, three Gal, one Neu5Ac: positioned on the C Gal.

- Three antennae, three Gal, two Neu5Ac: positioned on the A and C Gal.

- Four antennae, three Gal, one Neu5Ac: with Gal on A, B, and C, Neu5Ac
  is positioned on the C Gal.

- Four antennae, three Gal, two Neu5Ac: with Gal on A, B, and C, Neu5Ac
  residues are positioned on the A and C Gal.

- Four antennae, four Gal, two Neu5Ac: positioned on the A and C Gal.

Other distributions are left unchanged. Branch labels describe topology,
not linkage positions or drawing order. For three antennae, C is on the
Man with one antenna; A and B are interchangeable. For four antennae,
A/B and C/D are the two pairs. With three Gal, A/B is the pair with two
Gal, and C is the occupied antenna of the other pair.

Only concrete, topological structures are accepted. Use
[`glyrepr::remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.html)
explicitly before calling this function on intact or partial structures.
Missing values are preserved.

Only simple GlcNAc, Gal-GlcNAc, and Neu5Ac-Gal-GlcNAc antennae are
normalized. Non-N-glycans, other antenna counts, and structures with
extended or modified antennae, other sialic acids, or floating parts or
substituents are returned unchanged. Core fucose on the reducing-end
GlcNAc and a bisecting GlcNAc on the central Man are preserved.

## Examples

``` r
x <- glyparse::auto_parse(paste0(
  "Gal(b1-4)GlcNAc(b1-2)[GlcNAc(b1-6)]Man(a1-6)",
  "[GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
))
normalize_n_glycan(glyrepr::remove_linkages(x))
#> <glycan_structure[1]>
#> [1] Gal(??-?)GlcNAc(??-?)Man(??-?)[GlcNAc(??-?)[GlcNAc(??-?)]Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-
#> # Unique structures: 1
```
