# Calculate Product-Substrate Ratios

Calculate the ratio between the motif quantification of enzyme products
and substrates for glycosyltransferases and sulfotransferases. For
glycoproteomics data, ratios are calculated independently for each
glycosite.

## Usage

``` r
product_substrate_ratio(exp, enzymes)
```

## Arguments

- exp:

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object.

- enzymes:

  A character vector of glycosyltransferase or sulfotransferase names,
  or a list of their
  [`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)
  objects. Starter glycosyltransferases are not supported because their
  substrates are not glycans.

## Value

A plain
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with a `product_substrate_ratio` assay. For glycomics data, `rowData()`
contains `enzyme`. For glycoproteomics data, it contains `enzyme`,
`protein`, and `protein_site`. `colData()` and metadata are preserved
from `exp`.

## Details

Each rule contributes its product and acceptor motif once.
Quantifications are summed across all rules of an enzyme before
division, including when multiple rules contain the same motif. Rule
`rejects` and `requires` are ignored. Ratios with zero substrate
quantification are returned as `NA`.

Motifs are matched strictly for intact glycan structures and leniently
for partial or reduced structures, consistent with other glyenzy
motif-based functions.

## References

Bao B, Kellman BP, Chiang AWT, et al. (2021). Correcting for sparsity
and interdependence in glycomics by accounting for glycan biosynthesis.
*Nature Communications*, 12, 4988.
[doi:10.1038/s41467-021-25183-5](https://doi.org/10.1038/s41467-021-25183-5)

## Examples

``` r
exp <- glyexp::real_experiment2[seq_len(10), seq_len(3)]
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
product_substrate_ratio(exp, "ST6GAL1")
#> class: SummarizedExperiment 
#> dim: 1 3 
#> metadata(2): exp_type glycan_type
#> assays(1): product_substrate_ratio
#> rownames(1): ST6GAL1
#> rowData names(1): enzyme
#> colnames(3): S1 S2 S3
#> colData names(1): group
```
