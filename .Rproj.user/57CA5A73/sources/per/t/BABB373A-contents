#' Predicts new response, for an object of class \code{isgl}.
#'
#' @param isgl Object of class \code{isgl}.
#' @param newX New observations, as supplied to the isgl algorithm.
#' @return Vector of predictions.
#'
#' @seealso isgl
predict.isgl = function(isgl, newX){
  X = t(t(newX)- isgl$X.transform$X.means)
  X = t(t(X) / isgl$X.transform$X.scale)

  if(isgl$type == "linear"){
    y = X%*%isgl$beta + isgl$intercept
    return(y)
  }else if(isgl$type == "logit"){
    model_prediction = X%*%isgl$beta + isgl$intercept
    p = (1 + exp( -model_prediction ))^(-1)
    return(p)
  }
}
