# find DIANN files

find DIANN files

## Usage

``` r
get_files_DIANN(path)
```

## Arguments

- path:

  directory path to search for DIANN output files

## Value

list with report TSV and FASTA file paths.

## Examples

``` r
td <- tempdir()
writeLines("protein", file.path(td, "database.fasta"))
writeLines("report", file.path(td, "report.tsv"))
get_files_DIANN(td)
#> $reporttsv
#> [1] "report.tsv"
#> 
#> $fasta
#> [1] "database.fasta"
#> 
```
