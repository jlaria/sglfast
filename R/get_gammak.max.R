#' Computes an upper bound for the group lasso regularization penalty.
#'
#' Computes an upper bound of \eqn{\gamma_k} so that the solution is non zero.
#'
#' @param data A list with x and y
#' @param group.length A vector of group lengths
#' @param k A value between 1 and the number of groups
#' @param type String: \code{'linear'} or \code{'logit'}
#' @param lambda A vector of length 2, with the values of \eqn{\lambda_1} and \eqn{\lambda_2}
#' @param beta A vector, with an approximate \eqn{\beta}
#' @param beta0 The intercept
#' @return A number. The upper bound for \eqn{\gamma_k}.

get_gammak.max = function(data, group.length, k, type, lambda, beta = rep(0, ncol(data$x)),
                          beta0 = 0){
  if( type == "linear" ){
    group.start= c(0, cumsum(group.length)) + 1
    indices = group.start[k]:(group.start[k]+group.length[k] - 1)
    z = 1/nrow(data$x)*t(data$x[,indices])%*%(data$y - as.matrix(data$x[,-indices])%*%as.matrix(beta[-indices]) - beta0)
    S = 0
    for(i in 1:group.length[k]){
      S = S + max(abs(z[i]) - lambda[1], 0)^2
    }
    gammaj.max = sqrt(S)/lambda[2]
  }else if( type == "logit"){
    group.start= c(0, cumsum(group.length)) + 1
    indices = group.start[k]:(group.start[k]+group.length[k] - 1)
    d = rep(0, nrow(data$x))
    for(i in 1:length(d)){
      d[i] = - (1 + exp(-beta0 - data$x[i,-indices]%*%beta[-indices]))^(-1) + data$y[i]
    }
    z = 1/nrow(data$x)*t(data$x[,indices])%*%d
    S = 0
    for(i in 1:group.length[k]){
      S = S + max(abs(z[i]) - lambda[1], 0)^2
    }
    gammaj.max = sqrt(S)/lambda[2]
  }
  return(gammaj.max)
}
