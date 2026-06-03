# DIANN helper functions

DIANN helper functions

## Usage

``` r
read_DIANN_output(
  diann.path,
  fasta.file,
  nrPeptides = 2,
  q_value = 0.01,
  isUniprot = TRUE,
  rev = "REV_"
)
```

## Arguments

- diann.path:

  path to diann-output.tsv

- fasta.file:

  path to fasta file

- nrPeptides:

  peptide threshold

- q_value:

  q value threshold

- isUniprot:

  logical, is the fasta file from UniProt

- rev:

  prefix for reversed/decoy entries

## Value

This function now errors and points users to prolfquapp.

## Examples

``` r
try(read_DIANN_output("report.tsv", "database.fasta"), silent = TRUE)
```
