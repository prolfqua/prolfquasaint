# SAINTexpress contrast analysis facade

SAINTexpress contrast analysis facade

SAINTexpress contrast analysis facade

## Value

An R6 class generator.

## Details

Encapsulates the SAINTexpress pipeline so that it can be reached through
the same facade-registry dispatch as the prolfqua built-in methods (lm,
limma, firth, ...). Wraps
[`protein_2localSaint`](https://prolfqua.github.io/prolfquasaint/reference/saintExpress.md)
-\>
[`runSaint`](https://prolfqua.github.io/prolfquasaint/reference/saintExpress.md)
-\>
[`ContrastsSAINTexpress`](https://prolfqua.github.io/prolfquasaint/reference/ContrastsSAINTexpress.md),
with bait / control / gene / protein-length columns resolved from the
annotation data attached to the input `LFQData`.

Differs from LM-style facades in three ways:

- The contrast set is derived from the bait column (one bait vs.
  controls), so the `contrasts` argument is ignored.

- The output schema uses `Bait` / `log2_EFCs` / `SaintScore` / `BFDR`;
  `p.value` is `NA`.

- The facade carries SAINT input tables and the raw SAINTexpress result
  via `extra_artifacts()` for downstream report rendering.

## See also

Other modelling:
[`ContrastsSAINTexpress`](https://prolfqua.github.io/prolfquasaint/reference/ContrastsSAINTexpress.md)

## Super class

[`prolfqua::ContrastsInterface`](https://wolski.github.io/prolfqua/reference/ContrastsInterface.html)
-\> `ContrastsSAINTFacade`

## Public fields

- `contrast`:

  wrapped
  [`ContrastsSAINTexpress`](https://prolfqua.github.io/prolfquasaint/reference/ContrastsSAINTexpress.md)
  object

- `saint_input`:

  SAINTexpress input tables (`inter`/`prey`/`bait`)

- `saint_result`:

  raw SAINTexpress result list

- `protein_id`:

  hierarchy key the SAINT `Prey` column is mapped back to

- `.lfqdata`:

  stored reference to input LFQData

- `.row_annot`:

  stored reference to row annotation

## Methods

### Public methods

- [`ContrastsSAINTFacade$new()`](#method-ContrastsSAINTFacade-new)

- [`ContrastsSAINTFacade$get_contrasts()`](#method-ContrastsSAINTFacade-get_contrasts)

- [`ContrastsSAINTFacade$get_Plotter()`](#method-ContrastsSAINTFacade-get_Plotter)

- [`ContrastsSAINTFacade$to_wide()`](#method-ContrastsSAINTFacade-to_wide)

- [`ContrastsSAINTFacade$get_missing()`](#method-ContrastsSAINTFacade-get_missing)

- [`ContrastsSAINTFacade$get_ora()`](#method-ContrastsSAINTFacade-get_ora)

- [`ContrastsSAINTFacade$get_rank()`](#method-ContrastsSAINTFacade-get_rank)

- [`ContrastsSAINTFacade$extra_artifacts()`](#method-ContrastsSAINTFacade-extra_artifacts)

- [`ContrastsSAINTFacade$clone()`](#method-ContrastsSAINTFacade-clone)

Inherited methods

- [`prolfqua::ContrastsInterface$column_description()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-column_description)
- [`prolfqua::ContrastsInterface$contrast_summary_table()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-contrast_summary_table)
- [`prolfqua::ContrastsInterface$filter_significant()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-filter_significant)
- [`prolfqua::ContrastsInterface$get_config()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_config)
- [`prolfqua::ContrastsInterface$get_contrast_sides()`](https://prolfqua.github.io/prolfqua/html/ContrastsInterface.html#method-ContrastsInterface-get_contrast_sides)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ContrastsSAINTFacade$new(
      lfqdata,
      modelstr = NULL,
      contrasts = NULL,
      row_annot = NULL,
      spc = FALSE,
      engine = "r"
    )

#### Arguments

- `lfqdata`:

  aggregated
  [`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.html)
  object with bait / control columns in the annotation

- `modelstr`:

  ignored — SAINT derives contrasts from the bait column; accepted for
  facade uniformity

- `contrasts`:

  ignored for the same reason

- `row_annot`:

  data.frame with protein-level annotation (must contain the hierarchy
  key column). When `NULL`, falls back to `lfqdata$hierarchy_keys()`
  content from the long table.

- `spc`:

  if `TRUE` run SAINTexpress-spc; default `FALSE` (intensity model).

- `engine`:

  SAINTexpress engine; one of `"r"` (default) or whatever
  [`runSaint`](https://prolfqua.github.io/prolfquasaint/reference/saintExpress.md)
  accepts.

------------------------------------------------------------------------

### Method `get_contrasts()`

get contrast results from the wrapped
[`ContrastsSAINTexpress`](https://prolfqua.github.io/prolfquasaint/reference/ContrastsSAINTexpress.md);
adds a `facade` column matching the lm/limma facade convention.

#### Usage

    ContrastsSAINTFacade$get_contrasts(...)

#### Arguments

- `...`:

  passed to `ContrastsSAINTexpress$get_contrasts`

------------------------------------------------------------------------

### Method `get_Plotter()`

plotter for SAINT contrasts (delegates)

#### Usage

    ContrastsSAINTFacade$get_Plotter(...)

#### Arguments

- `...`:

  passed to `ContrastsSAINTexpress$get_Plotter`

------------------------------------------------------------------------

### Method `to_wide()`

wide format (delegates)

#### Usage

    ContrastsSAINTFacade$to_wide(...)

#### Arguments

- `...`:

  passed to `ContrastsSAINTexpress$to_wide`

------------------------------------------------------------------------

### Method `get_missing()`

protein × bait pairs absent from SAINT output

#### Usage

    ContrastsSAINTFacade$get_missing()

#### Returns

data.frame with hierarchy and bait columns

------------------------------------------------------------------------

### Method `get_ora()`

ORA input table (delegates)

#### Usage

    ContrastsSAINTFacade$get_ora(
      up = TRUE,
      FDR_threshold = 0.05,
      diff_threshold = 1
    )

#### Arguments

- `up`:

  if TRUE return bait-enriched prey

- `FDR_threshold`:

  BFDR threshold

- `diff_threshold`:

  log2 EFC threshold

------------------------------------------------------------------------

### Method `get_rank()`

rank table for enrichment tools (delegates; default score is `log2_EFCs`
because SAINTexpress has no p-value)

#### Usage

    ContrastsSAINTFacade$get_rank(score = "log2_EFCs")

#### Arguments

- `score`:

  column to use as rank score

------------------------------------------------------------------------

### Method `extra_artifacts()`

SAINT-specific artifacts to surface in downstream reports. Returns a
named list with the SAINT input tables and the raw result `list`.

#### Usage

    ContrastsSAINTFacade$extra_artifacts()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastsSAINTFacade$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
