#' DIANN helper functions
#' @param diann.path path to diann-output.tsv
#' @param fasta.file path to fasta file
#' @param nrPeptides peptide threshold
#' @param q_value q value threshold
#' @param isUniprot logical, is the fasta file from UniProt
#' @param rev prefix for reversed/decoy entries
#' @importFrom rlang .data
#' @export
#' @examples
#' \dontrun{
#' x <- prolfquapp::get_DIANN_files("inst/application/DIANN/2517219/")
#' xd <- read_DIANN_output(x$data, x$fasta)
#' debug(read_DIANN_output)
#' }
read_DIANN_output <- function(
  diann.path,
  fasta.file,
  nrPeptides = 2,
  q_value = 0.01,

  isUniprot = TRUE,
  rev = "REV_"
) {
  stop(
    "read_DIANN_output() moved out of prolfquasaint because DIANN parsing ",
    "belongs to the prolfquapp report facade. Use prolfquapp::preprocess_DIANN() ",
    "or prolfquapp::diann_read_output() instead.",
    call. = FALSE
  )
}


#' set factors and sample names columns
#' @param atable an AnalysisConfiguration object
#' @param msdata data.frame with annotation columns
#' @param repeated logical, look for subject/BioReplicate column
#' @param SAINT logical, use Bait_ instead of Group_ as factor name
#' @export
#' @examples
#'
#' dataset_csv_examples <- system.file("application/dataset_csv", package = "prolfquapp")
#' files <- dir(dataset_csv_examples, pattern = "*.csv")
#'
#' res <- list()
#' for (i in seq_along(files)) {
#'   res[[i]] <- readr::read_csv(file.path(dataset_csv_examples, files[i]),
#'     show_col_types = FALSE)
#' }
#'
#' for (i in seq_along(res)) {
#'   atable <- prolfqua::AnalysisConfiguration$new()
#'   atable$hierarchy[["protein_Id"]] <- c("Protein")
#'   atable$hierarchy[["peptide_Id"]] <- c("Peptide")
#'   tmp <- prolfquasaint::dataset_set_factors_deprecated(atable, res[[i]])
#'   cat(i, " : ", length(tmp$atable$factors), "factors : ",
#'     paste(tmp$atable$factors, collapse = "; "), "\n")
#' }
#'
dataset_set_factors_deprecated <- function(
  atable,
  msdata,
  repeated = TRUE,
  SAINT = FALSE
) {
  if (sum(grepl("^name", colnames(msdata), ignore.case = TRUE)) > 0) {
    atable$sample_name <- grep(
      "^name",
      colnames(msdata),
      value = TRUE,
      ignore.case = TRUE
    )
  }

  stopifnot(
    sum(grepl(
      "^channel|^Relative|^raw",
      colnames(msdata),
      ignore.case = TRUE
    )) >=
      1
  )
  fileName <- grep(
    "^channel|^Relative|^raw",
    colnames(msdata),
    value = TRUE,
    ignore.case = TRUE
  )[1]
  atable$file_name <- fileName

  stopifnot(
    sum(grepl(
      "^group|^bait|^Experiment",
      colnames(msdata),
      ignore.case = TRUE
    )) >=
      1
  )

  groupingVAR <- grep(
    "^group|^bait|^Experiment",
    colnames(msdata),
    value = TRUE,
    ignore.case = TRUE
  )
  if (any(grepl("^bait", groupingVAR, ignore.case = TRUE))) {
    groupingVAR <- grep("^bait", groupingVAR, value = TRUE, ignore.case = TRUE)[
      1
    ]
  } else {
    groupingVAR <- groupingVAR[1]
  }

  msdata[[groupingVAR]] <- gsub("[[:space:]]", "", msdata[[groupingVAR]])
  msdata[[groupingVAR]] <- gsub("[-\\+\\/\\*]", "_", msdata[[groupingVAR]])

  if (SAINT) {
    atable$factors[["Bait_"]] = groupingVAR
  } else {
    atable$factors[["Group_"]] = groupingVAR
  }

  atable$factor_depth <- 1

  if (
    sum(grepl(
      "^subject|^BioReplicate",
      colnames(msdata),
      ignore.case = TRUE
    )) ==
      1 &
      repeated
  ) {
    subvar <- grep(
      "^subject|^BioReplicate",
      colnames(msdata),
      value = TRUE,
      ignore.case = TRUE
    )
    atable$factors[["Subject_"]] = subvar

    fct <- dplyr::distinct(msdata[, c(atable$file_name, groupingVAR, subvar)])
    tmp <- data.frame(table(fct[, c(groupingVAR, subvar)]))
    if (all(tmp$Freq > 1)) {
      atable$factor_depth <- 2
    }
  }
  if (sum(grepl("^control", colnames(msdata), ignore.case = TRUE)) == 1) {
    atable$factors[["CONTROL"]] = grep(
      "^control",
      colnames(msdata),
      value = TRUE,
      ignore.case = TRUE
    )
  }
  return(list(atable = atable, msdata = msdata))
}
