#' Solves the sparse-group lasso problem
#'
#' Computes the solution of the SGL problem, defined as,
#' \deqn{\hat{\beta}(\lambda, \gamma)=\arg\min_{\beta\in B} \left\{ \hat{R}(\beta) + \lambda_2 \sum_{j=1}^{J} \gamma_j \|\beta^{(j)}\|_2 +  \lambda_1 \|\beta\|_1 \right\}.}
#'
#' @param data.train A list with matrix \code{x} and response vector \code{y}.
#' @param group.length A vector of group lengths. Columns of \code{x} are in the same order of groups.
#' @param lambdas A vector of length \eqn{J+2}. The first two components are those of \eqn{\lambda}, and the remaining components are the group weights \eqn{\gamma}.
#' @param type A string. Whether \code{'linear'} or \code{'logit'}
#' @param betas.old A vector, with an approximate initial solution \eqn{\beta}. Be careful!
#' @param simple Boolean, indicating whether we are solving the two-parameter version.
#' @return A list with the fit \code{beta} and \code{intercept}.
#'
solve_inner_problem <- function(data.train, group.length, lambdas, type = "linear",
                                betas.old = rep(0, ncol(data.train$x)), simple = F ){
  if( simple ){
    Sol = SGL(data = data.train,
              group.length = group.length,
              type = type,
              lambda1 = lambdas[1],
              lambda2 = lambdas[2])
  }else{
    Sol = SGL(data = data.train,
              group.length = group.length,
              type = type,
              lambda1 = lambdas[1],
              lambda2 = lambdas[2],
              groupW = lambdas[3:length(lambdas)],
              betas.old = betas.old)
  }
  return(Sol)
}
