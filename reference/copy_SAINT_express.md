# copy Markdown and runscript for DIANN diann-output.tsv

copy Markdown and runscript for DIANN diann-output.tsv

## Usage

``` r
copy_SAINT_express(workdir = getwd(), run_script = FALSE)
```

## Arguments

- workdir:

  directory where to copy file - default is current working directory.

- run_script:

  if TRUE, also copy the R run scripts

## Value

paths of the copied files; a file that cannot be copied gives a warning
and is left out.

## Examples

``` r
copy_SAINT_express(workdir = tempdir(), run_script = FALSE)
#> [1] "/tmp/RtmpHCqDCI/bibliography.bib"               
#> [2] "/tmp/RtmpHCqDCI/SaintExpressReportMsFragger.Rmd"
```
