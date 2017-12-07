#' Standarize the data for the sparse-group lasso
#'
#' @param data.train A list with matrix x and response y
#' @param data.validate A list with matrix x and response y
#' @return A list with the transformed data, the transformation,
#' and intercept.
standarize_inital_data = function(data.train, data.validate, type="linear"){

  X.transform = NULL
  intercept = 0
  X <- data.train$x
  means <- apply(X,2,mean)
  X <- t(t(X) - means)
  var <- apply(X,2,function(x)(sum(x^2)/(length(x)-1)))
  X <- t(t(X) / var)
  data.train$x <- X

  X.transform <- list(X.scale = var, X.means = means)

  X = data.validate$x
  X = t(t(X)-means)
  X = t(t(X) / var)
  data.validate$x = X

  if( type == "linear" ){
    intercept = mean(data.train$y)
    data.train$y = data.train$y - intercept
    data.validate$y = data.validate$y - intercept
  }

  return( list(data.train = data.train, data.validate = data.validate,
               X.transform = X.transform, intercept = intercept) )
}
