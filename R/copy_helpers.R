#' copy SAINTexpress doc file
#' @param workdir directory where to copy file - default is current working directory.
#' @keywords internal
#' @export
#' @examples
#' copy_SAINTe_doc(workdir = tempdir())
copy_SAINTe_doc <- function(workdir = getwd()){
  runscripts <- c("SaintExpress/SAINTexpress-manual.docx")
  .scriptCopyHelperVec(runscripts, workdir = workdir)
}

.run_markdown_with_params <-
  function(params,
           markdown_path,
           dest_path,
           dest_file_name,
           workdir = tempdir(),
           packagename = "prolfqua",
           format = "pdf") {
    res <- .scriptCopyHelperVec(markdown_path,
                                workdir = workdir,
                                packagename = packagename)
    dist_file_path <-
      file.path(dest_path, paste0(dest_file_name, ".", format))
    if (is.null(res)) {
      return(NULL)
    }
    rmarkdown::render(
      res[1],
      output_format = if (format == "pdf") {
        bookdown::pdf_document2()
      } else{
        bookdown::html_document2()
      },
      params = params,
      envir = new.env()
    )

    pdf_doc <- paste0(tools::file_path_sans_ext(res[1]), ".", format)
    message("XXXX--------------------------------------XXXX")
    if (pdf_doc != dist_file_path) {
      message("from " , pdf_doc, " to ", dist_file_path)
      file.copy(pdf_doc, dist_file_path, overwrite = TRUE)
    }
    return(dist_file_path)
  }
