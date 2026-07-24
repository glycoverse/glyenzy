# Get all enzymes

Return a named list of all built-in enzymes, or a character vector of
gene symbols if `return_str` is `TRUE`.

## Usage

``` r
db_enzymes(
  return_str = FALSE,
  include_starter_gt = TRUE,
  include_npre_gt = TRUE
)
```

## Arguments

- return_str:

  If `FALSE` (default), returns the enzyme list. Otherwise returns a
  character vector of gene symbols.

- include_starter_gt:

  If `TRUE` (default), includes starter GTs in the result. Starter GTs
  are enzymes that initiate glycosylation by introducing the first sugar
  residue onto a non-glycan substrate. For example, DPAGT1 is the
  starter GT for N-glycosylation.

- include_npre_gt:

  If `TRUE` (default), includes GTs involved in N-glycan precursor
  synthesis in the result. These GTs are responsible for building the
  N-glycan precursor before it is transferred to the target protein by
  OST.

## Value

A list of
[`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)s
or a character vector.

## Examples

``` r
db_enzymes(return_str = TRUE)
#>   [1] "ALG13"      "ALG14"      "ALG1"       "ALG2"       "ALG11"     
#>   [6] "ALG3"       "ALG9"       "ALG12"      "ALG6"       "ALG8"      
#>  [11] "ALG10"      "MOGS"       "GANAB"      "MAN1B1"     "MAN1A1"    
#>  [16] "MAN1A2"     "MAN1C1"     "MAN2A1"     "MAN2A2"     "MGAT1"     
#>  [21] "MGAT2"      "MGAT3"      "MGAT4A"     "MGAT4B"     "MGAT5"     
#>  [26] "MGAT5B"     "B3GNT2"     "B3GNT3"     "B3GNT4"     "B3GNT5"    
#>  [31] "B3GNT6"     "B3GNT7"     "B3GNT8"     "LFNG"       "MFNG"      
#>  [36] "RFNG"       "B3GALT1"    "B3GALT2"    "B3GALT4"    "B3GALT5"   
#>  [41] "B3GALT6"    "B4GALT1"    "B4GALT2"    "B4GALT3"    "B4GALT4"   
#>  [46] "B4GALT7"    "FUT1"       "FUT2"       "FUT3"       "FUT4"      
#>  [51] "FUT5"       "FUT6"       "FUT7"       "FUT8"       "FUT9"      
#>  [56] "ST6GAL1"    "ST6GAL2"    "ST3GAL1"    "ST3GAL2"    "ST3GAL3"   
#>  [61] "ST3GAL4"    "ST3GAL6"    "ST6GALNAC1" "ST6GALNAC2" "ST6GALNAC3"
#>  [66] "ST6GALNAC4" "ST6GALNAC5" "ST6GALNAC6" "ST8SIA1"    "ST8SIA2"   
#>  [71] "ST8SIA3"    "ST8SIA4"    "ST8SIA5"    "ST8SIA6"    "A4GNT"     
#>  [76] "C1GALT1"    "GCNT1"      "GCNT2"      "GCNT3"      "GCNT4"     
#>  [81] "B3GALNT1"   "B3GALNT2"   "B4GALNT1"   "B4GALNT2"   "B4GALNT3"  
#>  [86] "B4GALNT4"   "CSGALNACT1" "CSGALNACT2" "B3GAT1"     "B3GAT2"    
#>  [91] "B3GAT3"     "CHPF"       "CHPF2"      "CHSY1"      "CHSY3"     
#>  [96] "EXT1"       "EXT2"       "EXTL1"      "EXTL2"      "EXTL3"     
#> [101] "HAS1"       "HAS2"       "HAS3"       "B3GLCT"     "A4GALT"    
#> [106] "ABO"        "GXYLT1"     "GXYLT2"     "LARGE1"     "LARGE2"    
#> [111] "XXYLT1"     "POMGNT1"    "POMGNT2"    "DPAGT1"     "FUT10"     
#> [116] "FUT11"      "POFUT1"     "POFUT2"     "POGLUT1"    "POGLUT2"   
#> [121] "POGLUT3"    "POMT1"      "POMT2"      "TMTC1"      "TMTC2"     
#> [126] "TMTC3"      "TMTC4"      "GALNT1"     "GALNT2"     "GALNT3"    
#> [131] "GALNT4"     "GALNT5"     "GALNT6"     "GALNT7"     "GALNT8"    
#> [136] "GALNT9"     "GALNT10"    "GALNT11"    "GALNT12"    "GALNT13"   
#> [141] "GALNT14"    "GALNT15"    "GALNT16"    "GALNT17"    "GALNT18"   
#> [146] "GALNT19"    "CHST1"      "CHST2"      "CHST3"      "CHST4"     
#> [151] "CHST5"      "CHST6"      "CHST8"      "CHST9"      "CHST10"    
#> [156] "GAL3ST2"    "GAL3ST3"    "GAL3ST4"   
```
