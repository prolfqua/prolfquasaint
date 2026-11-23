#' copy SAINTexpress doc file
#' @param workdir directory where to copy file - default is current working directory.
#' @return Invisibly returns the copied file path, or `character(0)` if the
#'   companion package/manual is unavailable.
#' @keywords internal
#' @export
#' @examples
#' copy_SAINTe_doc(workdir = tempdir())
copy_SAINTe_doc <- function(workdir = getwd()) {
  if (!requireNamespace("saintexpressbin", quietly = TRUE)) {
    warning(
      "copy_SAINTe_doc requires the 'saintexpressbin' package, which is not installed."
    )
    return(invisible(character(0)))
  }
  src <- system.file(
    "manual",
    "SAINTexpress-manual.docx",
    package = "saintexpressbin"
  )
  if (!nzchar(src)) {
    warning("SAINTexpress-manual.docx not found in saintexpressbin.")
    return(invisible(character(0)))
  }
  dest <- file.path(workdir, basename(src))
  file.copy(src, dest, overwrite = TRUE)
  invisible(dest)
}


#' copy Markdown and runscript for DIANN diann-output.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @param run_script if TRUE, also copy the R run scripts
#' @return paths of the copied files; a file that cannot be copied gives a warning and is left out.
#' @export
#' @examples
#' copy_SAINT_express(workdir = tempdir(), run_script = FALSE)
copy_SAINT_express <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/bibliography.bib",
    "application/SE2/SaintExpressReportMsFragger.Rmd",
    if (run_script) {
      "application/SE2/CreateSaintExpress_Report.R"
    },
    if (run_script) {
      "application/SE2_DIANN/DIANN_SE.R"
    }
  )
  .copy_package_files(runscripts, workdir = workdir)
}

#' Copy files shipped in prolfquasaint into a working directory
#'
#' The report template is rendered in the working directory, because `rmarkdown::render()` writes next to its
#' input and the installed package is not writable. prolfquapp has the same helper, but prolfquasaint does not
#' import prolfquapp (prolfquapp loads prolfquasaint).
#' @param files paths relative to the installed package directory
#' @param workdir directory to copy the files to
#' @return paths of the copied files
#' @noRd
.copy_package_files <- function(files, workdir) {
  copied <- character()
  for (file in files) {
    src <- system.file(file, package = "prolfquasaint")
    dest <- file.path(workdir, basename(file))
    if (!nzchar(src)) {
      warning(sprintf("Could not copy %s: not found in prolfquasaint.", file), call. = FALSE)
    } else if (file.copy(src, dest, overwrite = TRUE)) {
      copied <- c(copied, dest)
    } else {
      warning(sprintf("Could not copy %s to %s.", src, dest), call. = FALSE)
    }
  }
  copied
}
