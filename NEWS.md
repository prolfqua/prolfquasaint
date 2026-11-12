# prolfquasaint 0.99.0

* The pure-R engine (`runSaint(engine = "r")`) comes from the package SAINTexpress, the new name of saintexpress since its 1.0.0 CRAN release.

* Installing prolfquasaint from GitHub (`pak::pak("prolfqua/prolfquasaint")`, CI) finds SAINTexpress in the `prolfqua/saintexpress` repository again; after the rename the install stopped with "Can't find package called SAINTexpress".

* `runSaint(engine = "r", optimizer = "nloptr")` stops with an error when nloptr is not installed, as saintexpress 0.99.4 does, instead of using `optim` without notice.

* The example datasets `data_SAINTe_output` and `data_IonstarProtein_subsetNorm` now ship with prolfquasaint (they were removed from prolfqua 1.9.0); load them with `prolfqua::prolfqua_data(<name>, package = "prolfquasaint")`.

* SAINT contrast tables now carry `estimate_type = "observed"`, like every other
  modelling backend. SAINTexpress imputes nothing, so the column was absent
  entirely, and a downstream consumer asking to exclude imputed estimates had no
  way to tell "nothing was imputed" from "provenance unknown".

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
