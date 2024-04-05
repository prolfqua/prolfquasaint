#' get apparams from bfabric executable
#' @export
apparams_Bfabric <- function(yml) {
  REPORTDATA <- list()

  # Application parameters
  REPORTDATA$spc <- FALSE
  REPORTDATA$FCthreshold <- if(!is.null(as.numeric( yml$application$parameters$`22|FCthreshold` ))){
    as.numeric( yml$application$parameters$`22|FCthreshold` ) } else { 2 }
  REPORTDATA$FDRthreshold <- if(!is.null(as.numeric( yml$application$parameters$`21|BFDRsignificance` ))){
    as.numeric(yml$application$parameters$`21|BFDRsignificance`) } else {0.1}
  REPORTDATA$Normalization <- yml$application$parameters$`11|Normalization`
  REPORTDATA$Transformation <- yml$application$parameters$`51|Transformation`
  nrPeptides <- yml$application$parameters$`61|nrPeptides`
  REPORTDATA$nrPeptides <- if(!is.null( nrPeptides )) { as.numeric(nrPeptides) } else { 2 }
  if(is.null( nrPeptides )) {warning("no prameter nrPeptides in yaml setting to 2")}
  return(REPORTDATA)
}

#' get params from bfabric executable
#' @export
get_params_Bfabric <- function(yml){
  BFABRIC <- list()
  BFABRIC$workunitID = yml$job_configuration$workunit_id
  BFABRIC$workunitURL = paste0("https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=",BFABRIC$workunitID,"&tab=details")
  BFABRIC$orderID = yml$job_configuration$order_id
  BFABRIC$inputID = purrr::map_chr(yml$job_configuration$input[[1]], "resource_id")
  BFABRIC$inputID = tail(BFABRIC$inputID,n = 1)
  BFABRIC$inputURL = purrr::map_chr(yml$job_configuration$input[[1]], "resource_url")
  BFABRIC$inputURL = tail(BFABRIC$inputURL, n = 1)

  BFABRIC$datasetID <- yml$application$parameters$`10|datasetId`
  return(BFABRIC)
}

#' normalize and then exponentiate data.
#' @export
normalize_exp <- function(lfqdataProt, normalization = c("vsn","robscale","none")) {
  normalization <- match.arg(normalization)
  if (normalization == "vsn") {
    lfqTrans <- prolfquapp::transform_lfqdata(lfqdataProt, method = "vsn")
    tr <- lfqTrans$get_Transformer()
    tr$intensity_array(exp, force = TRUE)
    tr$lfq$config$table$is_response_transformed <- FALSE
    lfqdataProt <- tr$lfq
  } else if (normalization == "robscale") {
    lfqTrans <- prolfquapp::transform_lfqdata(lfqdataProt, method = "robscale")
    tr <- lfqTrans$get_Transformer()
    tr$intensity_array(exp, force = TRUE)
    tr$lfq$config$table$is_response_transformed <- FALSE
    lfqdataProt <- tr$lfq
  } else {
    lfqdataProt <- lfqdataProt
  }
  return(lfqdataProt)
}


#' force data transformation
#' @export
transform_force <- function(lfqdataProt,  transformation = c("sqrt","log2","none")) {
  transformation <- match.arg(transformation)
  if (transformation == "sqrt") {
    tr <- lfqdataProt$get_Transformer()
    tr$intensity_array(sqrt, force = TRUE)
    tr$lfq$config$table$is_response_transformed <- FALSE
    lfqdataProt <- tr$lfq
  } else if (transformation == "log2") {
    tr <- lfqdataProt$get_Transformer()
    tr$intensity_array(log2, force = TRUE)
    tr$lfq$config$table$is_response_transformed <- FALSE
    lfqdataProt <- tr$lfq
  } else {

  }
  return(lfqdataProt)
}
