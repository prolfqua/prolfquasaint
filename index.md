# prolfquasaint

`prolfquasaint` is the [prolfqua](https://github.com/fgcz/prolfqua)
integration layer for SAINTexpress-based protein interaction scoring. It
converts prolfqua `LFQData` objects into SAINT input tables, runs
SAINTexpress, and exposes the results through `prolfqua` contrast
adapters.

📖 Documentation: <https://prolfqua.github.io/prolfquasaint/>

## Architecture

The SAINTexpress concerns are split across three packages:

- [`saintexpress`](https://github.com/prolfqua/saintexpress) — pure-R
  implementation of the spectral-count and intensity scoring engines.
- [`saintexpressbin`](https://github.com/prolfqua/saintexpressbin) —
  ships the native `SAINTexpress-spc` and `SAINTexpress-int` binaries
  plus a Docker fallback on macOS.
- `prolfquasaint` (this package) — prolfqua integration. Depends on
  `saintexpress` (`Imports`); `saintexpressbin` is `Suggests`. If
  `engine = "binary"` is requested without `saintexpressbin` installed,
  the call falls back to the R engine with a warning.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("prolfqua/saintexpress")
remotes::install_github("prolfqua/saintexpressbin")  # optional, for native engine
remotes::install_github("fgcz/prolfqua")
remotes::install_github("prolfqua/prolfquasaint")
```

`prolfquasaint`’s `Remotes:` field references `prolfqua/saintexpress`
and `prolfqua/saintexpressbin`, so
`pak::pkg_install("prolfqua/prolfquasaint")` should resolve them
transitively. If your dependency resolver does not walk nested
`Remotes:`, install `saintexpress` and `saintexpressbin` explicitly
first (see
[TODO/TODO_transitive_remotes.md](https://prolfqua.github.io/prolfquasaint/TODO/TODO_transitive_remotes.md)).

## Public API

The user-facing surface is stable:

- [`protein_2localSaint()`](https://prolfqua.github.io/prolfquasaint/reference/saintExpress.md)
  — convert a tidy long-format protein quantification table into SAINT
  `inter / prey / bait` data frames.
- `runSaint(si, ..., engine = c("binary", "r"), optimizer = c("base", "nloptr"))`
  — run SAINTexpress. `engine = "binary"` delegates to
  [`saintexpressbin::saintexpress_run()`](https://rdrr.io/pkg/saintexpressbin/man/saintexpress_run.html);
  `engine = "r"` delegates to
  [`saintexpress::run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.html).
- `ContrastsSAINTexpress`, `ContrastsSAINTFacade` — prolfqua R6 adapters
  that expose SAINT results through the prolfqua contrast interface.

## TIP49 Reference Comparison

`inst/test/` holds the TIP49 reference dataset and the comparison
harness that gates the R engine against the native SAINTexpress 3.6.3
binary. The `testthat` suite includes the cross-package comparison:

``` bash
make test
```

`inst/test/run_tip49_comparison.sh` regenerates the comparison reports
(`tip49_comparison_363_spc.html` and the matching diff). Use it after
changes to either engine to refresh the on-disk evidence.

## Vignettes

- `native-r-spc-comparison.qmd` — spectral-count engine: native binary
  vs R implementation on the TIP49 fixture.
- `native-r-int-comparison.qmd` — intensity engine: native binary vs R
  implementation on the TIP49 fixture.

Build with `make build-vignettes`.

For introductory walk-throughs with simulated data, see the vignettes
shipped by [`saintexpress`](https://github.com/prolfqua/saintexpress)
and [`saintexpressbin`](https://github.com/prolfqua/saintexpressbin).

## SAINTexpress Upstream

The original SAINTexpress 3.6.3 source archive was downloaded from the
`saint-apms` SourceForge project:

- SourceForge files page:
  <https://sourceforge.net/projects/saint-apms/files/>
- Archive: `SAINTexpress_v3.6.3__2018-03-09.tar.gz`
- SourceForge listing date: 2018-03-15

The extracted and cleaned C++ source lives in a separate repository;
this package depends only on prebuilt binaries shipped by
`saintexpressbin`.
