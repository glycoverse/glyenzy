# Get all enzymes

Return a named list of all built-in enzymes, or a character vector of
gene symbols if `return_str` is `TRUE`.

## Usage

``` r
all_enzymes(return_str = FALSE)
```

## Arguments

- return_str:

  If `FALSE` (default), returns the enzyme list. Otherwise returns a
  charactor vector of gene symbols.

## Value

A list of
[`enzyme()`](https://glycoverse.github.io/glyenzy/reference/enzyme.md)s
or a character vector.

## Examples

``` r
all_enzymes(return_str = TRUE)
#>   [1] "MOGS"       "GANAB"      "MAN1B1"     "MAN1A1"     "MAN1A2"    
#>   [6] "MAN1C1"     "MAN2A1"     "MAN2A2"     "MGAT1"      "MGAT2"     
#>  [11] "MGAT3"      "MGAT4A"     "MGAT4B"     "MGAT5"      "MGAT5B"    
#>  [16] "B3GNT2"     "B3GNT3"     "B3GNT4"     "B3GNT5"     "B3GNT6"    
#>  [21] "B3GNT7"     "B3GNT8"     "LFNG"       "MFNG"       "RFNG"      
#>  [26] "B3GALT1"    "B3GALT2"    "B3GALT4"    "B3GALT5"    "B3GALT6"   
#>  [31] "B4GALT1"    "B4GALT2"    "B4GALT3"    "B4GALT4"    "B4GALT5"   
#>  [36] "B4GALT6"    "B4GALT7"    "FUT1"       "FUT2"       "FUT3"      
#>  [41] "FUT4"       "FUT5"       "FUT6"       "FUT7"       "FUT8"      
#>  [46] "FUT9"       "ST6GAL1"    "ST6GAL2"    "ST3GAL1"    "ST3GAL2"   
#>  [51] "ST3GAL3"    "ST3GAL4"    "ST3GAL5"    "ST3GAL6"    "ST6GALNAC1"
#>  [56] "ST6GALNAC2" "ST6GALNAC3" "ST6GALNAC4" "ST6GALNAC5" "ST6GALNAC6"
#>  [61] "ST8SIA1"    "ST8SIA2"    "ST8SIA3"    "ST8SIA4"    "ST8SIA5"   
#>  [66] "ST8SIA6"    "A4GNT"      "C1GALT1"    "GCNT1"      "GCNT2"     
#>  [71] "GCNT3"      "GCNT4"      "B3GALNT1"   "B3GALNT2"   "B4GALNT1"  
#>  [76] "B4GALNT2"   "B4GALNT3"   "B4GALNT4"   "CSGALNACT1" "CSGALNACT2"
#>  [81] "B3GAT1"     "B3GAT2"     "B3GAT3"     "CHPF"       "CHPF2"     
#>  [86] "CHSY1"      "CHSY3"      "EXT1"       "EXT2"       "EXTL1"     
#>  [91] "EXTL2"      "EXTL3"      "HAS1"       "HAS2"       "HAS3"      
#>  [96] "B3GLCT"     "A4GALT"     "ABO"        "GXYLT1"     "GXYLT2"    
#> [101] "LARGE1"     "LARGE2"     "XXYLT1"     "POMGNT1"    "POMGNT2"   
```
