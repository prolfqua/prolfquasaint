#' copy SAINTexpress doc file
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#' @examples
#' copy_SAINTe_doc(workdir = tempdir())
copy_SAINTe_doc <- function(workdir = getwd()){
  runscripts <- c("SaintExpress/SAINTexpress-manual.docx")
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir)
}


#' copy Markdown and runscript for DIANN diann-output.tsv
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_SAINT_express <- function(workdir = getwd(), run_script = FALSE) {
  runscripts <- c(
    "application/bibliography.bib",
    "application/SE2/SaintExpressReportMsFragger.Rmd",
    if (run_script) {"application/SE2/CreateSaintExpress_Report.R"},
    if (run_script) {"application/SE2_DIANN/DIANN_SE.R"}
  )
  prolfqua:::.scriptCopyHelperVec(runscripts, workdir = workdir, packagename = "prolfquapp")
}

