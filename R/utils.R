#' get apparams from bfabric executable
#' @param yml parsed YAML configuration list from B-Fabric
#' @return list of SAINTexpress application parameters.
#' @export
#' @examples
#' yml <- list(application = list(parameters = list(
#'   `22|FCthreshold` = "2",
#'   `21|BFDRsignificance` = "0.1",
#'   `11|Normalization` = "none",
#'   `51|Transformation` = "none",
#'   `61|nrPeptides` = "2"
#' )))
#' apparams_Bfabric(yml)
apparams_Bfabric <- function(yml) {
  REPORTDATA <- list()

  # Application parameters
  REPORTDATA$spc <- FALSE
  REPORTDATA$FCthreshold <- if (
    !is.null(as.numeric(yml$application$parameters$`22|FCthreshold`))
  ) {
    as.numeric(yml$application$parameters$`22|FCthreshold`)
  } else {
    2
  }
  REPORTDATA$FDRthreshold <- if (
    !is.null(as.numeric(yml$application$parameters$`21|BFDRsignificance`))
  ) {
    as.numeric(yml$application$parameters$`21|BFDRsignificance`)
  } else {
    0.1
  }
  REPORTDATA$Normalization <- yml$application$parameters$`11|Normalization`
  REPORTDATA$Transformation <- yml$application$parameters$`51|Transformation`
  nrPeptides <- yml$application$parameters$`61|nrPeptides`
  REPORTDATA$nrPeptides <- if (!is.null(nrPeptides)) {
    as.numeric(nrPeptides)
  } else {
    2
  }
  if (is.null(nrPeptides)) {
    warning("no prameter nrPeptides in yaml setting to 2")
  }
  return(REPORTDATA)
}

#' get params from bfabric executable
#' @param yml parsed YAML configuration list from B-Fabric
#' @return list with B-Fabric workunit, order, and input identifiers.
#' @export
#' @examples
#' yml <- list(
#'   job_configuration = list(
#'     workunit_id = "1",
#'     order_id = "2",
#'     input = list(list(list(resource_id = "3", resource_url = "https://example.org")))
#'   ),
#'   application = list(parameters = list(`10|datasetId` = "4"))
#' )
#' get_params_Bfabric(yml)
get_params_Bfabric <- function(yml) {
  BFABRIC <- list()
  BFABRIC$workunitID <- yml$job_configuration$workunit_id
  BFABRIC$workunitURL <- paste0(
    "https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=",
    BFABRIC$workunitID,
    "&tab=details"
  )
  BFABRIC$orderID <- yml$job_configuration$order_id
  BFABRIC$inputID <- purrr::map_chr(
    yml$job_configuration$input[[1]],
    ~ as.character(.x$resource_id)
  )
  BFABRIC$inputID <- utils::tail(BFABRIC$inputID, n = 1)
  BFABRIC$inputURL <- purrr::map_chr(
    yml$job_configuration$input[[1]],
    "resource_url"
  )
  BFABRIC$inputURL <- utils::tail(BFABRIC$inputURL, n = 1)

  BFABRIC$datasetID <- yml$application$parameters$`10|datasetId`
  return(BFABRIC)
}

#' normalize and then exponentiate data.
#' @param lfqdataProt an LFQData object
#' @param normalization normalization method: "vsn", "robscale", or "none"
#' @return This function now errors and points users to prolfquapp.
#' @export
#' @examples
#' try(normalize_exp(NULL, normalization = "none"), silent = TRUE)
normalize_exp <- function(
  lfqdataProt,
  normalization = c("vsn", "robscale", "none")
) {
  normalization <- match.arg(normalization)
  stop(
    "normalize_exp() moved out of prolfquasaint because LFQData normalization ",
    "belongs to the prolfquapp report facade. Use prolfquapp::transform_lfqdata() ",
    "in prolfquapp before calling prolfquasaint::protein_2localSaint().",
    call. = FALSE
  )
}


#' force data transformation
#' @param lfqdataProt an LFQData object
#' @param transformation transformation method: "sqrt", "log2", or "none"
#' @return transformed LFQData object.
#' @export
#' @examples
#' obj <- list()
#' transform_force(obj, transformation = "none")
transform_force <- function(
  lfqdataProt,
  transformation = c("sqrt", "log2", "none")
) {
  transformation <- match.arg(transformation)
  if (transformation == "sqrt") {
    tr <- lfqdataProt$get_Transformer()
    tr$intensity_array(sqrt, force = TRUE)
    tr$lfq$config$is_response_transformed <- FALSE
    lfqdataProt <- tr$lfq
  } else if (transformation == "log2") {
    tr <- lfqdataProt$get_Transformer()
    tr$intensity_array(log2, force = TRUE)
    tr$lfq$config$is_response_transformed <- FALSE
    lfqdataProt <- tr$lfq
  } else {}
  return(lfqdataProt)
}

#' find DIANN files
#' @param path directory path to search for DIANN output files
#' @return list with report TSV and FASTA file paths.
#' @export
#' @examples
#' td <- tempdir()
#' writeLines("protein", file.path(td, "database.fasta"))
#' writeLines("report", file.path(td, "report.tsv"))
#' get_files_DIANN(td)
get_files_DIANN <- function(path) {
  diann.path <- grep(
    "report\\.tsv$|diann-output\\.tsv",
    dir(path = path, recursive = TRUE),
    value = TRUE
  )
  fasta.files <- grep(
    "*\\.fasta$|*\\.fas$",
    dir(path = path, recursive = TRUE),
    value = TRUE
  )

  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  if (length(fasta.files) == 0) {
    logger::log_error("No fasta file found!")
    stop()
  }

  return(list(reporttsv = diann.path, fasta = fasta.files))
}
