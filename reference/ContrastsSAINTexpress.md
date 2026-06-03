# Wrapper to results produced by SAINTexpress

Wrapper to results produced by SAINTexpress

Wrapper to results produced by SAINTexpress

## Value

An R6 class generator.

## Details

SAINTexpress writes

## See also

[`summary_ROPECA_median_p.scaled`](https://wolski.github.io/prolfqua/reference/summary_ROPECA_median_p.scaled.html)

Other modelling:
[`ContrastsSAINTFacade`](https://prolfqua.github.io/prolfquasaint/reference/ContrastsSAINTFacade.md)

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.html)
-\> `ContrastsSAINTexpress`

## Public fields

- `contrast_result`:

  data.frame with the contrast computation results

- `modelName`:

  model name

- `subject_id`:

  subject id defualt 'Prey'

## Methods

### Public methods

- [`ContrastsSAINTexpress$new()`](#method-ContrastsSAINTexpress-new)

- [`ContrastsSAINTexpress$get_contrast_sides()`](#method-ContrastsSAINTexpress-get_contrast_sides)

- [`ContrastsSAINTexpress$get_linfct()`](#method-ContrastsSAINTexpress-get_linfct)

- [`ContrastsSAINTexpress$get_contrasts()`](#method-ContrastsSAINTexpress-get_contrasts)

- [`ContrastsSAINTexpress$get_Plotter()`](#method-ContrastsSAINTexpress-get_Plotter)

- [`ContrastsSAINTexpress$to_wide()`](#method-ContrastsSAINTexpress-to_wide)

- [`ContrastsSAINTexpress$get_rank()`](#method-ContrastsSAINTexpress-get_rank)

- [`ContrastsSAINTexpress$get_ora()`](#method-ContrastsSAINTexpress-get_ora)

- [`ContrastsSAINTexpress$clone()`](#method-ContrastsSAINTexpress-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$extra_artifacts()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-extra_artifacts)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_config()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_missing()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_missing)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsSAINTexpress$new(
      contrastsdf,
      subject_id = "Prey",
      modelName = "ContrastSaint"
    )

#### Arguments

- `contrastsdf`:

  return value of
  [`runSaint`](https://prolfqua.github.io/prolfquasaint/reference/saintExpress.md)

- `subject_id`:

  default "Prey"

- `modelName`:

  name of model

------------------------------------------------------------------------

### Method `get_contrast_sides()`

show contrasts

#### Usage

    ContrastsSAINTexpress$get_contrast_sides()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_linfct()`

no available for SaintExpress

#### Usage

    ContrastsSAINTexpress$get_linfct()

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrasts

#### Usage

    ContrastsSAINTexpress$get_contrasts(all = FALSE)

#### Arguments

- `all`:

  should all columns be returned (default FALSE)

- `global`:

  use a global linear function (determined by get_linfct)

------------------------------------------------------------------------

### Method `get_Plotter()`

get
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.html)

#### Usage

    ContrastsSAINTexpress$get_Plotter(
      fc_threshold = 1,
      SaintScore = 0.75,
      bfdr_threshold = 0.1
    )

#### Arguments

- `fc_threshold`:

  fold change threshold to show

- `SaintScore`:

  SaintScore threshold to show in the heatmap.

- `bfdr_threshold`:

  BFDR threshold

#### Returns

[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.html)

------------------------------------------------------------------------

### Method `to_wide()`

convert to wide format

#### Usage

    ContrastsSAINTexpress$to_wide(columns = c("SaintScore", "BFDR"))

#### Arguments

- `columns`:

  value column default SaintScore, BFDR

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_rank()`

get signed rank-list input table for enrichment tools.

#### Usage

    ContrastsSAINTexpress$get_rank(score = "log2_EFCs")

#### Arguments

- `score`:

  column to use as rank score

#### Returns

data.frame with subject id, contrast, and score columns.

------------------------------------------------------------------------

### Method `get_ora()`

get SAINT ORA input table.

#### Usage

    ContrastsSAINTexpress$get_ora(
      up = TRUE,
      FDR_threshold = 0.05,
      diff_threshold = 1
    )

#### Arguments

- `up`:

  if TRUE return bait-enriched prey, otherwise negative effects

- `FDR_threshold`:

  BFDR threshold

- `diff_threshold`:

  log2 EFC threshold

#### Returns

filtered SAINT contrast data.frame for ORA input.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsSAINTexpress$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
seout <- prolfqua::prolfqua_data("data_SAINTe_output")
cse <- ContrastsSAINTexpress$new(seout$list)
stopifnot(dim(cse$to_wide()) == c(64,13))
cse$get_contrast_sides()
#>       contrast group_1 group_2
#> 1 b vs Control       b Control
#> 2 c vs Control       c Control
#> 3 d vs Control       d Control
#> 4 e vs Control       e Control
stopifnot(dim(cse$get_contrasts()) == c(236,7))
cse$get_linfct()
#> NULL
pl <- cse$get_Plotter()
stopifnot(c("gg", "ggplot") %in% class(pl$volcano()$FDR))
```
