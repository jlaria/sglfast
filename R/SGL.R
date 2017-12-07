#' Function to compute the SGL solution
#'
#' Originally written by Simon for R package SGL.

SGL <- function(data, group.length, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.1,
                gamma = 0.8, verbose = FALSE, step = 1, reset = 10,
                lambda1 = 0.2, lambda2 = 0.01, groupW = NULL, betas.old = rep(0, ncol(data$x))
){
  intercept = 0
  if(type == "linear"){
    Sol <- run_sgl_linear(data, group.length, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh,
                          min.frac = min.frac, lambda1 = lambda1, lambda2=lambda2,
                          gamma = gamma, verbose = verbose, step = step, reset = reset,
                          groupW = groupW, beta.naught = betas.old)
    Sol = list(beta = Sol$beta, intercept = intercept)
  }
  if(type == "logit"){
    Sol <- run_sgl_logit(data, group.length, thresh = thresh, inner.iter = maxit, outer.iter = maxit,
                         outer.thresh = thresh, min.frac = min.frac, lambda1 = lambda1, lambda2= lambda2,
                         gamma = gamma, verbose = verbose, step = step, reset = reset,
                         groupW = groupW, beta.naught = betas.old)
    Sol = list(beta = Sol$beta, intercept = Sol$intercept)
  }
  return(Sol)
}
