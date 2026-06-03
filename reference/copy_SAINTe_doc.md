# copy SAINTexpress doc file

copy SAINTexpress doc file

## Usage

``` r
copy_SAINTe_doc(workdir = getwd())
```

## Arguments

- workdir:

  directory where to copy file - default is current working directory.

## Value

Invisibly returns the copied file path, or \`character(0)\` if the
companion package/manual is unavailable.

## Examples

``` r
copy_SAINTe_doc(workdir = tempdir())
```
