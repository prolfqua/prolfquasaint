#' SAINTexpress example output
#'
#' Output of a SAINTexpress run, used by the examples of
#' \code{\link{ContrastsSAINTexpress}}. Moved here from prolfqua 1.9.0.
#'
#' @format A list of three data frames: `listFile` (path of the SAINTexpress
#'   `list.txt`), `list` (236 interaction scores, 17 columns) and `out` (the
#'   console output of the run, one line per row).
#' @family data
#' @docType data
#' @keywords internal
"data_SAINTe_output"

#' IonStar protein-level example data
#'
#' Median-polish protein abundances of a subset of the IonStar E. coli spike-in
#' dilution series (PRIDE PXD003881), normalized. Used by the examples of
#' \code{\link{protein_2localSaint}}. Moved here from prolfqua 1.9.0.
#'
#' @format A list with `data`, a data frame with 3260 rows (columns
#'   `protein_Id`, `sampleName`, `dilution.`, `run_Id`, `raw.file`, `isotope`,
#'   `medpolish`), and `config`, its \code{\link[prolfqua]{AnalysisConfiguration}}.
#' @source PRIDE PXD003881
#' @family data
#' @docType data
#' @keywords internal
"data_IonstarProtein_subsetNorm"
