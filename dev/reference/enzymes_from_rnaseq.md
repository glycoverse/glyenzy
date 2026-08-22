# Filter enzymes using RNA-seq expression

`enzymes_from_rnaseq()` selects built-in enzymes whose genes are
expressed at or above a TPM threshold. When `mat` contains multiple
columns, the mean TPM across columns is used for each gene.

## Usage

``` r
enzymes_from_rnaseq(mat, threshold = 1)
```

## Arguments

- mat:

  A numeric matrix of TPM values with gene symbols as unique row names
  and samples or replicates as columns.

- threshold:

  A non-negative numeric scalar giving the minimum mean TPM for a gene
  to be considered expressed. The default follows Huang et al. (2021):
  genes with TPM below 1 are treated as not expressed or rarely
  expressed, whereas genes with TPM at least 1 are treated as expressed.

## Value

A named list of
[`enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/enzyme.md)
objects. Genes not represented in the built-in enzyme database are
ignored.

## References

Huang, Y.-F., Aoki, K., Akase, S., et al. (2021). Global mapping of
glycosylation pathways in human-derived cells. *Developmental Cell*, 56,
1195-1209.
[doi:10.1016/j.devcel.2021.02.023](https://doi.org/10.1016/j.devcel.2021.02.023)
.

## Examples

``` r
tpm <- matrix(
  c(0.2, 1, 3),
  ncol = 1,
  dimnames = list(c("FUT8", "ST3GAL3", "MOGS"), "sample")
)
enzymes_from_rnaseq(tpm)
#> $MOGS
#> 
#> ── Enzyme: MOGS ────────────────────────────────────────────────────────────────
#> ℹ Type: "GH" (Glycoside hydrolase)
#> ℹ Species: "human"
#> ℹ Glycan type: "N"
#> 
#> ── Rules (1) ──
#> 
#> → Rule 1: whole alignment
#> Acceptor:
#> "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#> Product:
#> "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#> 
#> $ST3GAL3
#> 
#> ── Enzyme: ST3GAL3 ─────────────────────────────────────────────────────────────
#> ℹ Type: "GT" (Glycosyltransferase)
#> ℹ Species: "human"
#> ℹ Glycan type: "lipid" and "free"
#> 
#> ── Rules (4) ──
#> 
#> → Rule 1: terminal alignment
#> Acceptor: "Gal(b1-3)GlcNAc(b1-"
#> Product: "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
#> Rejects:
#> "Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-"
#> → Rule 2: terminal alignment
#> Acceptor: "Gal(b1-4)GlcNAc(b1-"
#> Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
#> Rejects:
#> "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-"
#> → Rule 3: core alignment
#> Acceptor: "Gal(b1-3)GalNAc(a1-"
#> Product: "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
#> → Rule 4: terminal alignment
#> Acceptor: "Gal(b1-4)Glc(b1-"
#> Product: "Neu5Ac(a2-3)Gal(b1-4)Glc(b1-"
#> 
```
