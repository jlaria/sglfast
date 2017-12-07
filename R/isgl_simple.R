#' Two-parameter version of the iterative sparse-group lasso
#'
#' Computes the solution of the sparse-group lasso problem, defined as,
#' \deqn{\hat{\beta}(\lambda, \gamma)=\arg\min_{\beta\in B} \left\{ \hat{R}(\beta) + \lambda_2 \sum_{j=1}^{J} \gamma_j \|\beta^{(j)}\|_2 +  \lambda_1 \|\beta\|_1 \right\},}
#' where \eqn{\gamma_j} is the square root of j-th group's size.
#' Regularization parameter \eqn{\lambda} is automatically selected, so that the validation error is minimized.
#'
#' @param data.train A list with the matrix of observations \code{x}, and the response \code{y}, to train the model.
#' @param data.validate A list with the matrix of observations \code{x}, and the response \code{y}, to validate the model.
#' @param index A vector with group indices.
#' @param group.length A vector with group length, alternatively to \code{index}.
#' @param type String: \code{'linear'} or \code{'logit'}.
#' @param standardize Whether to stardardize \code{data.train} and \code{data.validate}. Recommended if \code{type} is \code{'logit'}.
#' @return An object of class \code{isgl}.

isgl_simple = function( data.train, data.validate, index = NULL, group.length = NULL, type = "linear",
                        standardize = F){

  # We tranform the initial data
  if (standardize){
    temp = standarize_inital_data(data.train, data.validate, type)
    data.train = temp$data.train
    data.validate =  temp$data.validate
    X.transform = temp$X.transform
    intercept = temp$intercept
  }
  else{
    X.transform = list(X.means = rep(0, ncol(data.train$x)), X.scale = rep(1, ncol(data.train$x)))
    intercept = 0
  }

  # We sort and compute group lengths
  if( is.null(group.length) ){
    if( is.null(index) ){
      write("Error 1: You must provide valid group indices or lengths")
      return(1)
    }
    temp = index2group.length(index)
    group.length = temp$group.length
    ord = temp$ord
    data.train$x = data.train$x[, ord]
    data.validate$x = data.validate$x[, ord]
    unord = match(1:length(ord),ord)
  }else{
    unord = 1:ncol(data.train$x)
  }

  # Compute initial lambdas
  lambda.max = c(0, 0)
  gamma = rep(0, length(group.length))
  lambda.max[1] = get_lambda1.max(data.train, group.length, type = type)

  lambda.init = c(lambda.max[1]*0.1, 1)
  for (i in 1:length(gamma)) {
    gamma[i] = get_gammak.max(data.train, group.length, i, type, lambda.init)
  }

  lambda.max[2] = max(gamma/sqrt(group.length))

  # Start the iterative search
  nparams = length(lambda.init)
  best_lambdas <- lambda.max*0.1
  num_solves <- 0
  max_solves <- 5000

  # Compute initial model
  model_params = solve_inner_problem(data.train, group.length, best_lambdas, type, simple = T)
  best_cost = get_validation_cost(data.validate$x, data.validate$y, model_params, type)
  best_beta = model_params
  num_solves = num_solves+1

  # Set initial coordinate
  coord <- 1
  fixed = 0

  # Main loop
  while ((num_solves < max_solves) &&(fixed < nparams) ) {
    old_lambdas <- best_lambdas

    # No need to compute lambda1,2 max?

    direction = 1
    curr_lambdas = best_lambdas
    stepsize = c( best_lambdas[coord]*0.1,
                  best_lambdas[coord]*0.1 )
    while (direction < 3) {

      if(stepsize[direction] > 0){
        step = runif(1, stepsize[direction]*0.1, stepsize[direction])*(3-2*direction)
        curr_lambdas[coord] = best_lambdas[coord] + step
        model_params <- solve_inner_problem(data.train, group.length, curr_lambdas, type, simple = T)
        num_solves <- num_solves + 1
        cost <- get_validation_cost( data.validate$x, data.validate$y, model_params, type)
        if(best_cost - cost > 0.00001){
          best_cost = cost
          best_lambdas <- curr_lambdas
          best_beta = model_params

          stepsize[direction] = c( stepsize[1]*2,
                                   min(stepsize[2]*2, best_lambdas[coord])
          )[direction]
        }
        else{
          direction = direction + 1
        }
      }else{
        direction = direction + 1
      }
    }

    if(old_lambdas[coord] == best_lambdas[coord]){
      fixed = fixed + 1
    }else{
      fixed = 0
    }
    coord <- coord%%(nparams) +1
  }

  solution = list(best_lambdas=best_lambdas, num_solves=num_solves, best_cost=best_cost,
                  beta = best_beta$beta[unord],
                  intercept = best_beta$intercept + intercept,
                  X.transform = X.transform,
                  type = type)
  class(solution) = "isgl"
  return( solution )
}
