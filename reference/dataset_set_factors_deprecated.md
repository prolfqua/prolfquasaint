# set factors and sample names columns

set factors and sample names columns

## Usage

``` r
dataset_set_factors_deprecated(atable, msdata, repeated = TRUE, SAINT = FALSE)
```

## Arguments

- atable:

  an AnalysisConfiguration object

- msdata:

  data.frame with annotation columns

- repeated:

  logical, look for subject/BioReplicate column

- SAINT:

  logical, use Bait\_ instead of Group\_ as factor name

## Value

list with updated \`atable\` and \`msdata\`.

## Examples

``` r
dataset_csv_examples <- system.file("application/dataset_csv", package = "prolfquapp")
files <- dir(dataset_csv_examples, pattern = "*.csv")

res <- list()
for (i in seq_along(files)) {
  res[[i]] <- readr::read_csv(file.path(dataset_csv_examples, files[i]),
    show_col_types = FALSE)
}

for (i in seq_along(res)) {
  atable <- prolfqua::AnalysisConfiguration$new()
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["peptide_Id"]] <- c("Peptide")
  tmp <- prolfquasaint::dataset_set_factors_deprecated(atable, res[[i]])
  cat(i, " : ", length(tmp$atable$factors), "factors : ",
    paste(tmp$atable$factors, collapse = "; "), "\n")
}
```
