# ContrastsSAINTFacade -----

#' SAINTexpress contrast analysis facade
#'
#' Encapsulates the SAINTexpress pipeline so that it can be reached
#' through the same facade-registry dispatch as the prolfqua built-in
#' methods (lm, limma, firth, ...). Wraps
#' \code{\link{protein_2localSaint}} -> \code{\link{runSaint}} ->
#' \code{\link{ContrastsSAINTexpress}}, with bait / control / gene /
#' protein-length columns resolved from the annotation data attached to
#' the input \code{LFQData}.
#'
#' Differs from LM-style facades in three ways:
#' \itemize{
#'   \item The contrast set is derived from the bait column (one bait
#'         vs. controls), so the \code{contrasts} argument is ignored.
#'   \item The output schema uses \code{Bait} / \code{log2_EFCs} /
#'         \code{SaintScore} / \code{BFDR}; \code{p.value} is \code{NA}.
#'   \item The facade carries SAINT input tables and the raw
#'         SAINTexpress result via \code{extra_artifacts()} for
#'         downstream report rendering.
#' }
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
ContrastsSAINTFacade <- R6::R6Class(
  "ContrastsSAINTFacade",
  inherit = prolfqua::ContrastsInterface,
  public = list(
    #' @field contrast wrapped \code{\link{ContrastsSAINTexpress}} object
    contrast = NULL,
    #' @field saint_input SAINTexpress input tables
    #'   (\code{inter}/\code{prey}/\code{bait})
    saint_input = NULL,
    #' @field saint_result raw SAINTexpress result list
    saint_result = NULL,
    #' @field protein_id hierarchy key the SAINT \code{Prey} column is
    #'   mapped back to
    protein_id = character(),
    #' @field .lfqdata stored reference to input LFQData
    .lfqdata = NULL,
    #' @field .row_annot stored reference to row annotation
    .row_annot = NULL,
    #' @description
    #' initialize
    #' @param lfqdata aggregated \code{\link[prolfqua]{LFQData}} object
    #'   with bait / control columns in the annotation
    #' @param modelstr ignored — SAINT derives contrasts from the bait
    #'   column; accepted for facade uniformity
    #' @param contrasts ignored for the same reason
    #' @param row_annot data.frame with protein-level annotation (must
    #'   contain the hierarchy key column). When \code{NULL}, falls back
    #'   to \code{lfqdata$hierarchy_keys()} content from the long table.
    #' @param spc if \code{TRUE} run SAINTexpress-spc; default
    #'   \code{FALSE} (intensity model).
    #' @param engine SAINTexpress engine; one of \code{"r"} (default)
    #'   or whatever \code{\link{runSaint}} accepts.
    initialize = function(
      lfqdata,
      modelstr = NULL,
      contrasts = NULL,
      row_annot = NULL,
      spc = FALSE,
      engine = "r"
    ) {
      subject_id <- lfqdata$subject_id()
      hierarchy_keys <- lfqdata$hierarchy_keys()
      if (!identical(subject_id, hierarchy_keys)) {
        stop(
          "ContrastsSAINTFacade requires aggregated LFQData. ",
          "`lfqdata$subject_id()` must equal `lfqdata$hierarchy_keys()`.",
          call. = FALSE
        )
      }
      if (!is.null(modelstr) && length(modelstr) > 0 && nzchar(modelstr)) {
        warning(
          "ContrastsSAINTFacade ignores `modelstr`: contrasts are derived ",
          "from the bait column.",
          call. = FALSE
        )
      }
      if (!is.null(contrasts) && length(contrasts) > 0) {
        warning(
          "ContrastsSAINTFacade ignores `contrasts`: bait vs. control ",
          "contrasts are derived from the annotation.",
          call. = FALSE
        )
      }
      self$.lfqdata <- lfqdata

      if (is.null(row_annot)) {
        protein_id <- lfqdata$relevant_hierarchy_keys()[[1]]
        row_annot <- dplyr::distinct(
          lfqdata$data_long()[, protein_id, drop = FALSE]
        )
      }
      self$.row_annot <- row_annot

      prepared <- .saint_prepare_data(lfqdata, row_annot)
      saint_input <- protein_2localSaint(
        prepared$data,
        quantcolumn = prepared$response,
        proteinID = prepared$protein_id,
        geneNames = ".saint_gene_names",
        proteinLength = ".saint_protein_length",
        IP_name = prepared$sample_col,
        baitCol = prepared$bait_col,
        CorTCol = prepared$control_col
      )
      saint_result <- runSaint(
        saint_input,
        spc = spc,
        engine = engine,
        CLEANUP = TRUE
      )
      # Remap Prey -> project hierarchy key so downstream joins on the
      # hierarchy key work without renaming.
      saint_result$list[[prepared$protein_id]] <- saint_result$list$Prey

      self$saint_input <- saint_input
      self$saint_result <- saint_result
      self$protein_id <- prepared$protein_id
      self$contrast <- ContrastsSAINTexpress$new(
        saint_result$list,
        subject_id = prepared$protein_id
      )
      self$config <- self$contrast$get_config()
    },
    #' @description get contrast results from the wrapped
    #'   \code{\link{ContrastsSAINTexpress}}; adds a \code{facade}
    #'   column matching the lm/limma facade convention.
    #' @param ... passed to \code{ContrastsSAINTexpress$get_contrasts}
    get_contrasts = function(...) {
      res <- self$contrast$get_contrasts(...)
      if (!("facade" %in% colnames(res))) {
        res <- dplyr::mutate(res, facade = "saint", .before = 1)
      }
      res
    },
    #' @description plotter for SAINT contrasts (delegates)
    #' @param ... passed to \code{ContrastsSAINTexpress$get_Plotter}
    get_Plotter = function(...) self$contrast$get_Plotter(...),
    #' @description wide format (delegates)
    #' @param ... passed to \code{ContrastsSAINTexpress$to_wide}
    to_wide = function(...) self$contrast$to_wide(...),
    #' @description protein × bait pairs absent from SAINT output
    #' @return data.frame with hierarchy and bait columns
    get_missing = function() {
      protein_id <- self$protein_id
      all_subjects <- dplyr::distinct(
        self$.lfqdata$data_long()[, protein_id, drop = FALSE]
      )
      contrast_sides <- self$contrast$get_contrast_sides()
      baits <- unique(contrast_sides$group_1)
      expected <- tidyr::crossing(all_subjects, Bait = baits)
      estimated <- dplyr::distinct(
        self$contrast$contrast_result[, c(protein_id, "Bait"), drop = FALSE]
      )
      dplyr::anti_join(expected, estimated, by = c(protein_id, "Bait"))
    },
    #' @description ORA input table (delegates)
    #' @param up if TRUE return bait-enriched prey
    #' @param FDR_threshold BFDR threshold
    #' @param diff_threshold log2 EFC threshold
    get_ora = function(up = TRUE, FDR_threshold = 0.05, diff_threshold = 1) {
      self$contrast$get_ora(
        up = up,
        FDR_threshold = FDR_threshold,
        diff_threshold = diff_threshold
      )
    },
    #' @description rank table for enrichment tools (delegates; default
    #'   score is \code{log2_EFCs} because SAINTexpress has no p-value)
    #' @param score column to use as rank score
    get_rank = function(score = "log2_EFCs") {
      self$contrast$get_rank(score = score)
    },
    #' @description SAINT-specific artifacts to surface in
    #'   downstream reports. Returns a named list with the SAINT input
    #'   tables and the raw result \code{list}.
    extra_artifacts = function() {
      list(
        saint_inter = self$saint_input$inter,
        saint_prey = self$saint_input$prey,
        saint_bait = self$saint_input$bait,
        saint_list = self$saint_result$list
      )
    }
  )
)
