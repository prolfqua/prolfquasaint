# TODO: Refactor SAINTexpress Support into Minimal Packages

## Goal

Split the current `prolfquasaint` responsibilities into two small standalone R
packages plus the existing prolfqua integration package:

- `saintexpressbin`: minimal package for shipping SAINTexpress native
  executables and resolving their paths.
- `saintexpress`: minimal package for the pure-R implementation and SAINT input
  table validation.
- `prolfquasaint`: prolfqua ecosystem integration, reports, fixtures, and
  comparison vignettes.

Design pressure: keep the two new leaf packages dependency-light and
independent of each other. `prolfquasaint` is the integration layer and
depends on both.

## Resolved Decisions

- **Package independence:** `saintexpress` is pure R. `saintexpressbin` ships
  binaries. Neither depends on the other. `prolfquasaint` depends on both.
- **Binary distribution:** `saintexpressbin` ships prebuilt binaries only.
  Source code for SAINTexpress lives in a separate GitHub repository and is
  not built at install time.
- **Native execution boundary:** Docker support and binary file execution live
  in `saintexpressbin`. `prolfquasaint` prepares SAINT input data frames, but
  `saintexpressbin` is responsible for writing `inter/prey/bait` files, running
  the native executable or Docker container, reading `list.txt`, and cleaning up
  native-run artifacts.
- **TIP49 fixtures:** live in `prolfquasaint` only. The leaf packages use
  minimal synthetic inputs for their own unit tests.
- **Comparison vignettes:** stay in `prolfquasaint` (the only package that
  depends on both leaves and holds the fixtures).
- **Compatibility API:** `prolfquasaint::runSaint()`,
  `prolfquasaint::protein_2localSaint()`, `ContrastsSAINTexpress`, and
  `ContrastsSAINTFacade` remain the stable public API indefinitely. Internally
  they delegate to `saintexpress::run_saint()` and `saintexpressbin`. No
  deprecation.

## Current State

`prolfquasaint` currently carries three concerns:

1. SAINTexpress executable and binary management.
2. Pure R implementations of `SAINTexpress-spc` and `SAINTexpress-int`.
3. prolfqua-facing integration, including `LFQData` conversion, contrast
   facades, reports, and comparison vignettes.

Relevant current files:

- `R/saint_spc.R`: spectral-count R implementation.
- `R/saint_int.R`: intensity R implementation.
- `R/saint_prepare.R`: SAINT input preparation from prolfqua objects.
- `R/tidyMS_SaintExpress.R`: SAINT input conversion, native execution helpers,
  and `runSaint()`-related support.
- `R/ContrastSaintExpress.R`: `ContrastsSAINTexpress` prolfqua contrast wrapper.
- `R/ContrastsSAINTFacade.R`: prolfqua modelling facade.
- `inst/SaintExpress/bin/`: packaged native executables.
- `inst/test/`: TIP49 fixtures and native-vs-reference comparison outputs.
- `vignettes/native-r-spc-comparison.qmd` and
  `vignettes/native-r-int-comparison.qmd`: comparison vignettes.

## Package Responsibilities

### `saintexpressbin`

Purpose: ship and locate SAINTexpress native executables.

Responsibilities:

- Ship prebuilt `SAINTexpress-spc` and `SAINTexpress-int` binaries under
  `inst/bin/<platform>/`.
- Resolve the platform-specific executable path.
- Provide a small runner that writes SAINT input files, shells out to the native
  binary or Docker, reads `list.txt`, and optionally cleans up generated files.
- Own Docker support for platforms where native binaries are unavailable or when
  Docker execution is explicitly requested.

Public API (proposed):

- `saintexpress_executable(type = c("spc", "int"))`
- `saintexpress_available(type = c("spc", "int"))`
- `saintexpress_run(si, type = c("spc", "int"), workdir = getwd(),
  cleanup = TRUE, use_docker = NULL)`

Dependencies: keep minimal, but use normal R package dependencies when they
avoid unnecessary reimplementation. No `prolfqua` and no `R6`.

Notes:

- Source for SAINTexpress lives in a separate GitHub repository. This package
  ships prebuilt binaries only; it does not build from source at install.
- Use explicit platform directories:
  - `inst/bin/Darwin/SAINTexpress-spc`
  - `inst/bin/Darwin/SAINTexpress-int`
  - `inst/bin/Linux64/SAINTexpress-spc`
  - `inst/bin/Linux64/SAINTexpress-int`
  - `inst/bin/Windows64/SAINTexpress-spc.exe`
  - `inst/bin/Windows64/SAINTexpress-int.exe`
- Ensure Unix binaries are executable in the package source and after install.
- Tests use injectable `sysname` for platform path resolution and minimal
  synthetic inputs (no TIP49 fixtures here).
- Verify the SAINTexpress redistribution/license situation before publishing
  this package publicly. Internal/private use is fine in the interim.

### `saintexpress`

Purpose: pure-R implementation of SAINTexpress algorithms and input table
validation.

Responsibilities:

- Validate SAINT input tables (`inter`, `prey`, `bait`).
- Run R implementations of `spc` and `int`.
- Return SAINTexpress-compatible `list.txt`-style data frames.

Public API (proposed):

- `run_saint(si, mode = c("spc", "int"), ...)` (R engine only)
- `validate_saint_input(si)`

Dependencies: keep modest, but do not force a pure-base rewrite of the existing
implementation. Use normal R package dependencies, including `dplyr` or optional
optimizer packages, if they materially reduce complexity. **No dependency on
`saintexpressbin`** — native dispatch is handled by `prolfquasaint` via
`saintexpressbin`, not here.

Move candidates from `prolfquasaint`:

- `R/saint_int.R`
- `R/saint_spc.R`
- shared parsing/model helpers used by both engines
- the R-engine side of `runSaint()`

Tests use minimal synthetic inputs (no TIP49 fixtures here).

### `prolfquasaint`

Purpose: prolfqua ecosystem integration layer.

Responsibilities:

- Convert `LFQData` and tidy prolfqua objects to SAINT input.
- Provide `ContrastsSAINTexpress` and `ContrastsSAINTFacade`.
- Depend on both `saintexpress` (R engine) and `saintexpressbin` (native
  binaries). This is the only package that bridges the two.
- Keep `runSaint()` and `protein_2localSaint()` as the stable public API,
  internally delegating to `saintexpress::run_saint()` for the R engine and
  `saintexpressbin::saintexpress_run()` for the native/Docker engine.
- Preserve the current `runSaint()` argument contract unless a later plan
  explicitly changes it: `filedir`, `spc`, `CLEANUP`, `use_docker`, `engine`,
  and `optimizer`.
- Prepare SAINT input data frames in `prolfquasaint`; pass those data frames to
  `saintexpressbin` for file writing and native execution.
- Own the TIP49 fixtures under `inst/test/`.
- Own the comparison vignettes (native binary vs R implementation,
  prolfqua/prolfquapp integration examples).
- Own the TIP49 cross-package integration tests that gate
  `saintexpress::run_saint()` against the reference SAINTexpress outputs.
- Keep report templates and B-Fabric/prolfquapp helpers unless they belong in
  `prolfquapp`.

Dependencies: as today plus `saintexpress` and `saintexpressbin` via
`Imports`/`Remotes`.

## Migration Plan

### 1. Create `saintexpressbin`

- Scaffold a minimal R package in a new repo.
- Move `inst/SaintExpress/bin/*` and the manual.
- Move binary path resolution, Docker support, file writing, native execution,
  output reading, and cleanup helpers; strip `prolfqua` dependencies.
- Add tests for platform path resolution with injectable `sysname`.
- Add tests for writing synthetic SAINT input files and reading native-style
  `list.txt` output without requiring TIP49 fixtures.

### 2. Create `saintexpress`

- Scaffold a minimal R package in a new repo.
- Move the pure-R SAINT implementation and table validation.
- Provide `run_saint()` as the standalone public interface (R engine only).
- Add tests using minimal synthetic inputs.
- Keep dependencies modest, but do not rewrite working implementation code only
  to avoid a dependency.

### 3. Refactor `prolfquasaint`

- Add `saintexpress` and `saintexpressbin` to `Imports` (with `Remotes:`
  entries until the packages are published).
- Replace internal R algorithm calls with `saintexpress::run_saint()`.
- Replace native binary lookup/run, Docker, file writing, output reading, and
  cleanup code with `saintexpressbin`.
- Keep prolfqua-facing names stable:
  - `protein_2localSaint()`
  - `runSaint()` (dispatches to R or native via the two new packages)
  - `ContrastsSAINTexpress`
  - `ContrastsSAINTFacade`
- Keep TIP49 fixtures in `inst/test/`.
- Keep TIP49 comparison tests in `prolfquasaint`, and make them the
  cross-package gate for R implementation correctness against SAINTexpress
  reference outputs.
- Update comparison vignettes to call both new packages explicitly.

### 4. Clean Up Dependencies

- `saintexpressbin`: minimal dependencies; no `prolfqua` and no `R6`.
- `saintexpress`: modest dependencies allowed when useful; no `prolfqua` and no
  dependency on `saintexpressbin`.
- `prolfquasaint`: existing deps plus `saintexpress` and `saintexpressbin`.
- Remove algorithm and binary-management code from `prolfquasaint` after tests
  pass.

### 5. Ecosystem integration

- Add `saintexpressbin` and `saintexpress` to the root `Makefile` `installs`
  target ahead of `prolfquasaint`.
- Run `make sync-r-makefiles` after scaffolding so both new packages use the
  canonical Makefile.
- Audit reverse dependencies (`prolfquapp`, `prophosqua`,
  `prolfquappPTMreaders`, `ptm-pipeline` templates) for direct calls to
  `runSaint()` / `protein_2localSaint()`. Since the public API is preserved,
  no downstream changes should be required — verify rather than assume.

### 6. Verify

Install and test in dependency order:

```bash
make install      # in saintexpressbin
make install      # in saintexpress
make install      # in prolfquasaint
make test         # in each package
make check-fast   # in each package
```

Then rebuild the comparison vignettes in `prolfquasaint`.

## Naming Note

Prefer `saintexpress` and `saintexpressbin`. They match the executable/manual
terminology and are clearer than `saintxpress`. Verify no name collisions on
CRAN/Bioconductor/GitHub before publishing.

## Remaining Open Items

- Public release vs. private hosting of `saintexpressbin` pending SAINTexpress
  license review.
- Exact GitHub org/owner for the two new repos and for the upstream
  SAINTexpress source repo referenced by `saintexpressbin`.
