test_that("prolfquasaint does not import prolfquapp", {
  description_path <- system.file("DESCRIPTION", package = "prolfquasaint")
  if (!nzchar(description_path)) {
    description_path <- test_path("..", "..", "DESCRIPTION")
  }
  description <- read.dcf(description_path)
  imports <- description[1, "Imports"]
  expect_no_match(imports, "\\bprolfquapp\\b")
})

test_that("moved prolfquapp-coupled helpers fail clearly", {
  expect_error(
    read_DIANN_output("report.tsv", "database.fasta"),
    "moved out of prolfquasaint"
  )
  expect_error(
    normalize_exp(NULL, normalization = "none"),
    "moved out of prolfquasaint"
  )
})
