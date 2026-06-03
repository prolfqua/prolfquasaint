# get params from bfabric executable

get params from bfabric executable

## Usage

``` r
get_params_Bfabric(yml)
```

## Arguments

- yml:

  parsed YAML configuration list from B-Fabric

## Value

list with B-Fabric workunit, order, and input identifiers.

## Examples

``` r
yml <- list(
  job_configuration = list(
    workunit_id = "1",
    order_id = "2",
    input = list(list(list(resource_id = "3", resource_url = "https://example.org")))
  ),
  application = list(parameters = list(`10|datasetId` = "4"))
)
get_params_Bfabric(yml)
#> $workunitID
#> [1] "1"
#> 
#> $workunitURL
#> [1] "https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=1&tab=details"
#> 
#> $orderID
#> [1] "2"
#> 
#> $inputID
#> [1] "3"
#> 
#> $inputURL
#> [1] "https://example.org"
#> 
#> $datasetID
#> [1] "4"
#> 
```
