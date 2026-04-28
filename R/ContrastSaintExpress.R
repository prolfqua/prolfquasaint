# ContrastsSAINTexpress ----
#' Wrapper to results produced by SAINTexpress
#'
#' SAINTexpress writes
#'
#'
#'
#' @export
#' @family modelling
#'
#' @examples
#' seout <- prolfqua::prolfqua_data("data_SAINTe_output")
#' cse <- ContrastsSAINTexpress$new(seout$list)
#' stopifnot(dim(cse$to_wide()) == c(64,13))
#' cse$get_contrast_sides()
#' stopifnot(dim(cse$get_contrasts()) == c(236,7))
#' cse$get_linfct()
#' pl <- cse$get_Plotter()
#' stopifnot(c("gg", "ggplot") %in% class(pl$volcano()$BFDR))
ContrastsSAINTexpress <- R6::R6Class(
  "ContrastsSAINTexpress",
  inherit = prolfqua::ContrastsInterface, # or remove it?
  public = list(
    #' @field contrast_result data.frame with the contrast computation results
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_id subject id defualt 'Prey'
    subject_id = character(),
    #' @description
    #' initialize
    #' @param contrastsdf return value of \code{\link{runSaint}}
    #' @param subject_id default "Prey"
    #' @param modelName name of model
    initialize = function(contrastsdf,
                          subject_id = "Prey",
                          modelName = "ContrastSaint"
    ){


      self$contrast_result = contrastsdf
      self$subject_id = subject_id
      self$modelName = modelName

      if ( "AvgIntensity" %in% colnames(contrastsdf)) {
        self$contrast_result <- contrastsdf |> dplyr::mutate(log2_EFCs = log2(FoldChange),
                                                      avgAbd = log2(AvgIntensity),
                                                      modelName = modelName)

      }else{
        self$contrast_result <- contrastsdf |> dplyr::mutate(log2_EFCs = log2(FoldChange),
                                                      avgAbd = log2(AvgSpec) ,
                                                      modelName = modelName)

      }

    },
    #' @description
    #' show contrasts
    #' @return data.frame
    get_contrast_sides = function(){
      dd <- self$contrast_result
      baits <- unique(dd$Bait)
      tt <- data.frame(contrast = paste0(baits, " vs Control"),
                       group_1  = baits,
                       group_2 = "Control")
      return(tt)
    },
    #' @description
    #' no available for SaintExpress
    #'
    get_linfct = function(){
      NULL
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link[prolfqua]{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      res <- self$contrast_result |> dplyr::select(
        tidyselect::all_of(c(self$subject_id,
                 "modelName",
                 "Bait",
                 "avgAbd",
                 "log2_EFCs",
                 "SaintScore",
                 "BFDR"
        )))
      res
    },
    #' @description get \code{\link[prolfqua]{ContrastsPlotter}}
    #' @param fc_threshold fold change threshold to show
    #' @param SaintScore SaintScore threshold to show in the heatmap.
    #' @param bfdr_threshold BFDR threshold
    #' @return \code{\link[prolfqua]{ContrastsPlotter}}
    get_Plotter = function(fc_threshold = 1, SaintScore = 0.75, bfdr_threshold = 0.1){
      res <- prolfqua::ContrastsPlotter$new(
        self$contrast_result,
        subject_id = self$subject_id,
        fcthresh = fc_threshold,
        volcano = list(list(score = "BFDR", thresh = bfdr_threshold)),
        histogram = list(list(score = "BFDR", xlim = c(0,1,0.05)), list(score = "SaintScore", xlim = c(0,1,0.05))),
        score = list(list(score = "SaintScore", thresh = SaintScore )),
        modelName = "modelName",
        diff = "log2_EFCs",
        contrast = "Bait")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default SaintScore, BFDR
    #' @return data.frame
    to_wide = function(columns = c("SaintScore", "BFDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- prolfqua::pivot_model_contrasts_to_wide(contrast_minimal,
                                                     subject_id = self$subject_id,
                                                     columns = c("log2_EFCs", columns),
                                                     contrast = 'Bait')
      return(contrasts_wide)
    }
  ))
