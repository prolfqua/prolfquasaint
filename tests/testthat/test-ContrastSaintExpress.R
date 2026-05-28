test_that("ContrastsSAINTexpress provides rank and ORA input tables", {
  saint_list <- data.frame(
    Bait = c("B1", "B1", "B1", "B2"),
    Prey = c("P1", "P2", "P3", "P4"),
    PreyGene = c("G1", "G2", "G3", "G4"),
    AvgSpec = c(10, 20, 30, 40),
    FoldChange = c(4, 0.25, 3, 8),
    SaintScore = c(0.9, 0.8, 0.7, 0.95),
    BFDR = c(0.01, 0.01, 0.2, 0.03)
  )

  contrasts <- ContrastsSAINTexpress$new(saint_list)

  plotter <- contrasts$get_Plotter()
  volcano <- plotter$volcano()
  expect_named(volcano, "FDR")
  expect_s3_class(volcano$FDR, "ggplot")
  expect_equal(volcano$FDR$labels$y, "-log10(BFDR)")

  rank_table <- contrasts$get_rank()
  expect_named(rank_table, c("Prey", "contrast", "score"))
  expect_equal(rank_table$score, log2(saint_list$FoldChange))

  ora_up <- contrasts$get_ora(
    up = TRUE,
    FDR_threshold = 0.05,
    diff_threshold = 1
  )
  expect_equal(ora_up$Prey, c("P1", "P4"))
  expect_equal(ora_up$contrast, c("B1", "B2"))

  ora_down <- contrasts$get_ora(
    up = FALSE,
    FDR_threshold = 0.05,
    diff_threshold = 1
  )
  expect_equal(ora_down$Prey, "P2")
  expect_equal(ora_down$contrast, "B1")
})
