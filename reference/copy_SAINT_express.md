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

return value from \`prolfqua::script_copy_helper_vec()\`.

## Examples

``` r
copy_SAINT_express(workdir = tempdir(), run_script = FALSE)
#> copy /home/runner/work/_temp/Library/prolfquasaint/application/bibliography.bib to /tmp/Rtmpfoz8EO/bibliography.bib
#> copy /home/runner/work/_temp/Library/prolfquasaint/application/SE2/SaintExpressReportMsFragger.Rmd to /tmp/Rtmpfoz8EO/SaintExpressReportMsFragger.Rmd
#> your working directory now should contain: 2 new files:
#> [1] "/tmp/Rtmpfoz8EO/bibliography.bib"               
#> [2] "/tmp/Rtmpfoz8EO/SaintExpressReportMsFragger.Rmd"
```
