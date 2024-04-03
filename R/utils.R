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
