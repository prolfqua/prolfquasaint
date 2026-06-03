# normalize and then exponentiate data.

normalize and then exponentiate data.

## Usage

``` r
normalize_exp(lfqdataProt, normalization = c("vsn", "robscale", "none"))
```

## Arguments

- lfqdataProt:

  an LFQData object

- normalization:

  normalization method: "vsn", "robscale", or "none"

## Value

This function now errors and points users to prolfquapp.

## Examples

``` r
try(normalize_exp(NULL, normalization = "none"), silent = TRUE)
```
