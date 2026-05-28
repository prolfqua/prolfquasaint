#' copy SAINTexpress doc file
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#' @examples
#' copy_SAINTe_doc(workdir = tempdir())
copy_SAINTe_doc <- function(workdir = getwd()){
  if (!requireNamespace("saintexpressbin", quietly = TRUE)) {
    warning("copy_SAINTe_doc requires the 'saintexpressbin' package, which is not installed.")
    return(invisible(character(0)))
  }
  src <- system.file("manual", "SAINTexpress-manual.docx", package = "saintexpressbin")
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
#' @export
#'
copy_SAINT_express <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/bibliography.bib",
    "application/SE2/SaintExpressReportMsFragger.Rmd",
    if (run_script) {"application/SE2/CreateSaintExpress_Report.R"},
    if (run_script) {"application/SE2_DIANN/DIANN_SE.R"}
  )
  prolfqua::script_copy_helper_vec(runscripts, workdir = workdir, packagename = "prolfquasaint")
}

