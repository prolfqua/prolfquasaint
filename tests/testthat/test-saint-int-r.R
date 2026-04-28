test_that("R intensity engine matches SAINTexpress-int TIP49 reference", {
  fixture_dir <- system.file(
    "test/saintexpress-363-tip49-reference-int",
    package = "prolfquasaint"
  )
  skip_if(fixture_dir == "", "SAINTexpress intensity fixture is not installed")

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

  r_out <- runSaint(si, spc = FALSE, engine = "r", optimizer = "base")$list
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
    expect_lte(max(delta, na.rm = TRUE), 0.005, label = column)
  }

  odds_delta <- abs(comparison$OddsScore_reference - comparison$OddsScore_r)
  expect_lte(max(odds_delta, na.rm = TRUE), 0.1)

  bfdr_delta <- abs(comparison$BFDR_reference - comparison$BFDR_r)
  expect_lte(max(bfdr_delta, na.rm = TRUE), 0.02)

  fold_change_rel_delta <- abs(
    comparison$FoldChange_reference - comparison$FoldChange_r
  ) /
    pmax(1, abs(comparison$FoldChange_reference))
  expect_lte(max(fold_change_rel_delta, na.rm = TRUE), 0.07)
})

test_that("R intensity engine preserves SAINTexpress output structure", {
  fixture_dir <- system.file(
    "test/saintexpress-363-tip49-reference-int",
    package = "prolfquasaint"
  )
  skip_if(fixture_dir == "", "SAINTexpress intensity fixture is not installed")

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

  r_out <- runSaint(si, spc = FALSE, engine = "r")$list

  self <- r_out[r_out$Bait == "ACTR5" & r_out$Prey == "ACTR5", ]
  expect_equal(nrow(self), 1)
  expect_equal(self$Intensity, ".")

  expect_true(any(grepl("\\.", r_out$ctrlIntensity, fixed = FALSE)))
  expect_true(all(r_out$BFDR >= 0 & r_out$BFDR <= 1))
  expect_equal(min(r_out$BFDR[r_out$AvgP >= 0.999]), 0)
})
