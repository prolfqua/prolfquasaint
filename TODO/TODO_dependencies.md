# TODO: Fix prolfquasaint Dependency Direction

## Aim

This document records the dependency-cycle issue found while testing the
`prolfquapp` Docker image build with vignette installation enabled. The goal is
to make `prolfquasaint` installable as a dependency of `prolfquapp`, without
`prolfquasaint` importing or requiring `prolfquapp` at install time.

## Current Blocker

`prolfquapp` imports `prolfquasaint`.

The published GitHub `prolfqua/prolfquasaint` `main` branch still imports
`prolfquapp`. This creates a circular dependency:

```text
prolfquapp -> prolfquasaint -> prolfquapp
```

During the `prolfquapp` Docker build this fails in the dependency solver before
the package can even reach the vignette build step:

```text
github::prolfqua/prolfquasaint: Can't install dependency prolfquapp
prolfquapp: Can't find package called prolfqua, prolfquapp
```

The local `prolfquasaint/DESCRIPTION` already removes `prolfquapp` from
`Imports`, but that fix is not yet available to Docker while `prolfquapp`
installs `github::prolfqua/prolfquasaint`.

## Required Dependency Direction

The intended package direction is:

```text
prolfqua <- prolfquasaint <- prolfquapp
```

Rules:

- `prolfquasaint` may import `prolfqua`.
- `prolfquapp` may import `prolfquasaint`.
- `prolfquasaint` must not import, depend on, or require `prolfquapp`.
- Any code that needs `prolfquapp` helpers belongs in `prolfquapp`, not in
  `prolfquasaint`.
- Examples or documentation may mention `prolfquapp`, but package load,
  install, tests, and runtime APIs must not require it.

## Refactor Tasks

- Keep `prolfquapp` removed from `prolfquasaint/DESCRIPTION` `Imports`.
- Check `DESCRIPTION`, `NAMESPACE`, `R/`, `tests/`, `vignettes/`, and
  `inst/` for remaining hard dependencies on `prolfquapp`.
- Move any remaining `prolfquapp`-coupled data parsing or LFQ normalization
  helpers out of `prolfquasaint`.
- If old helper names must remain temporarily, keep them as clear erroring
  compatibility stubs only when this is an intentional API decision.
- Ensure `ContrastsSAINTexpress`, `ContrastsSAINTFacade`, `runSaint()`, and
  `protein_2localSaint()` work without `prolfquapp` installed or loaded.
- Keep the integration direction explicit in tests:
  `prolfquapp` tests may exercise `prolfquasaint`; `prolfquasaint` tests must
  not require `prolfquapp`.

## Tests and Acceptance Checks

Run these before considering the dependency refactor complete:

```bash
Rscript -e "devtools::test()"
R CMD build .
R CMD INSTALL prolfquasaint_*.tar.gz
```

Add or keep tests that assert:

- `DESCRIPTION` `Imports` does not contain `prolfquapp`.
- `library(prolfquasaint)` succeeds in a clean library where `prolfquapp` is
  not installed.
- `prolfquasaint` unit tests pass without loading `prolfquapp`.
- Public SAINT facade APIs remain available after the refactor.

After publishing or otherwise making the fixed `prolfquasaint` ref available,
rerun the `prolfquapp` Docker build:

```bash
docker build --progress=plain -t prolfquapp-vignettes-test .
```

The `prolfquapp` Docker build should then advance past `prolfquasaint`
dependency resolution and reach the vignette-build checks for:

```text
doc/Grp2Analysis_V2_R6.Rmd
doc/DiffExpQC_R6.Rmd
```

## Notes for prolfquapp Docker

The `prolfquapp` Dockerfile can install `prolfquasaint` explicitly only after
the dependency cycle is fixed upstream. Until then, adding
`github::prolfqua/prolfquasaint` to `pak::pkg_install()` fails because the
published package metadata still points back to `prolfquapp`.

Do not work around this in `prolfquapp` by skipping dependency resolution or
copying files manually. The root cause is the package dependency direction in
`prolfquasaint`.
