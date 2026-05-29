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

## Examples

``` r
if (FALSE) { # \dontrun{
x <- prolfquapp::get_DIANN_files("inst/application/DIANN/2517219/")
xd <- read_DIANN_output(x$data, x$fasta)
debug(read_DIANN_output)
} # }
```
