# get apparams from bfabric executable

get apparams from bfabric executable

## Usage

``` r
apparams_Bfabric(yml)
```

## Arguments

- yml:

  parsed YAML configuration list from B-Fabric

## Value

list of SAINTexpress application parameters.

## Examples

``` r
yml <- list(application = list(parameters = list(
  `22|FCthreshold` = "2",
  `21|BFDRsignificance` = "0.1",
  `11|Normalization` = "none",
  `51|Transformation` = "none",
  `61|nrPeptides` = "2"
)))
apparams_Bfabric(yml)
#> $spc
#> [1] FALSE
#> 
#> $FCthreshold
#> [1] 2
#> 
#> $FDRthreshold
#> [1] 0.1
#> 
#> $Normalization
#> [1] "none"
#> 
#> $Transformation
#> [1] "none"
#> 
#> $nrPeptides
#> [1] 2
#> 
```
