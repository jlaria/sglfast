#' Evaluates the risk function
#'
#' @param X The model matrix
#' @param y The response vector
#' @param model_params The sgl fit.
#' @param type String: \code{'linear'} or \code{'logit'}.
#' @return The sum of errors
get_validation_cost <- function(X, y, model_params, type = "linear"){
  if( type == 'linear' ){
    return(test_error_linear(X, y, model_params))
  }else if(type=='logit'){
    return(test_error_logit(X, y, model_params))
  }
}
