#' Computes the maximum weight of the lasso penalty
#'
#' Given the \code{data} and \eqn{\beta}, this function computes an upper bound for \eqn{\lambda_1} such that \eqn{\beta(\lambda_1)\ne 0}.
#'
#' @param data A list with x and y
#' @param group.length A vector with group lengths
#' @param beta A vector, the previous solution \eqn{\beta}
#' @param beta0 A number, the intercept
#' @param type A string. Type of loss function: \code{'logit'} or \code{'linear'}
#' @return A number. The maximum \eqn{\lambda_1}
get_lambda1.max = function(data, group.length, beta = rep(0, ncol(data$x)), beta0 = 0,  type){
  if (type == "linear"){
    d = rep(0, ncol(data$x))
    for(i in 1:length(d)){
      d[i] = abs(1/nrow(data$x) * data$x[,i]%*%(data$y - data$x[, -i]%*%beta[-i] - beta0) )
    }
  }else if(type == "logit"){
    if(any(beta!=0) || beta0!=0){
      d = rep(0, ncol(data$x))
      for(j in 1:length(d)){
        D = rep(0, nrow(data$x))
        for(i in 1:length(D)){
          D[i] = - (1 + exp(-beta0 - data$x[i,-j]%*%beta[-j]))^(-1) + data$y[i]
        }
        d[j] = abs( 1/nrow(data$x)*data$x[,j]%*%D )
      }
    }else{
      d = rep(0, ncol(data$x))
      for(j in 1:length(d)){
        D = - 0.5 + data$y
        d[j] = abs( 1/nrow(data$x)*data$x[,j]%*%D )
      }
    }
  }
  return(max(d))
}
