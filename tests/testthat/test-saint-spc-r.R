test_that("R spectral-count engine matches SAINTexpress-spc TIP49 reference", {
  fixture_dir <- system.file(
    "test/saintexpress-363-tip49-reference-spc",
    package = "prolfquasaint"
  )
  skip_if(
    fixture_dir == "",
    "SAINTexpress spectral-count fixture is not installed"
  )

  si <- list(
    inter = read.delim(
      file.path(fixture_dir, "inter.dat"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    prey = read.delim(
      file.path(fixture_dir, "prey.dat"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    bait = read.delim(
      file.path(fixture_dir, "bait.dat"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  r_out <- runSaint(si, spc = TRUE, engine = "r", optimizer = "base")$list
  reference <- read.delim(
    file.path(fixture_dir, "list.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_named(r_out, names(reference))
  expect_equal(nrow(r_out), nrow(reference))

  key <- c("Bait", "Prey", "PreyGene")
  comparison <- merge(
    reference,
    r_out,
    by = key,
    suffixes = c("_reference", "_r")
  )
  expect_equal(nrow(comparison), nrow(reference))

  probability_cols <- c("AvgP", "MaxP", "TopoAvgP", "TopoMaxP", "SaintScore")
  for (column in probability_cols) {
    delta <- abs(
      comparison[[paste0(column, "_reference")]] -
        comparison[[paste0(column, "_r")]]
    )
    expect_lte(max(delta, na.rm = TRUE), 0.011, label = column)
  }

  odds_delta <- abs(
    comparison$logOddsScore_reference - comparison$logOddsScore_r
  )
  expect_lte(max(odds_delta, na.rm = TRUE), 0.011)

  bfdr_delta <- abs(comparison$BFDR_reference - comparison$BFDR_r)
  expect_lte(max(bfdr_delta, na.rm = TRUE), 0.011)

  expect_equal(comparison$FoldChange_reference, comparison$FoldChange_r)
})

test_that("R spectral-count engine preserves SAINTexpress output structure", {
  fixture_dir <- system.file(
    "test/saintexpress-363-tip49-reference-spc",
    package = "prolfquasaint"
  )
  skip_if(
    fixture_dir == "",
    "SAINTexpress spectral-count fixture is not installed"
  )

  si <- list(
    inter = read.delim(
      file.path(fixture_dir, "inter.dat"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    prey = read.delim(
      file.path(fixture_dir, "prey.dat"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    bait = read.delim(
      file.path(fixture_dir, "bait.dat"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  r_out <- runSaint(si, spc = TRUE, engine = "r")$list

  self <- r_out[r_out$Bait == "ACTR5" & r_out$Prey == "ACTR5", ]
  expect_equal(nrow(self), 1)
  expect_equal(self$Spec, "0")

  expect_true(any(grepl("\\|", r_out$ctrlCounts)))
  expect_true(all(r_out$BFDR >= 0 & r_out$BFDR <= 1))
  expect_equal(min(r_out$BFDR[r_out$AvgP >= 0.99]), 0)
})
