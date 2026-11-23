.local_dir <- function(env = parent.frame()) {
  path <- tempfile("copy_helpers_")
  dir.create(path)
  do.call(on.exit, list(substitute(unlink(path, recursive = TRUE), list(path = path)), add = TRUE), envir = env)
  path
}

test_that("copy_SAINT_express copies the report template and its bibliography", {
  path <- .local_dir()
  copied <- copy_SAINT_express(workdir = path, run_script = FALSE)
  expect_setequal(basename(copied), c("bibliography.bib", "SaintExpressReportMsFragger.Rmd"))
  expect_true(all(file.exists(copied)))
})

test_that("a file missing from the package warns and is left out", {
  path <- .local_dir()
  expect_warning(copied <- .copy_package_files("application/no-such-file.Rmd", path), "not found in prolfquasaint")
  expect_length(copied, 0)
})
