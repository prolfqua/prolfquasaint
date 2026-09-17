# prolfquasaint 0.99.0

* Development installs and CI now follow the current `saintexpress` repository instead of requesting a pinned historical commit.

* `get_rank()` now accepts `score = NULL` to mean "use the effect size this
  backend reports", as the `ContrastsInterface` contract defines it. Passing
  `NULL` previously selected no column at all and produced a rank table with no
  score, which broke callers that let the backend choose its own rank.

# prolfquasaint 0.1.6

* Require the SAINTexpress control-variance fix so container installations cannot reuse an older
  `saintexpress` revision that fails on constant control profiles.

# prolfquasaint 0.1.5

* Declared `cyclocomp` as a development dependency so package linting works in a fresh checkout.
* Declared the non-CRAN dependencies `prolfqua` and `saintexpressbin` in `Remotes` so the package
  installs from a fresh clone (previously only `saintexpress` was listed).
* Prepared package metadata, documentation, vignettes, and CI for Bioconductor-style checks.
