run_sgl_logit <- function(data, group.length,  thresh = 0.0001, lambda1 = 0.2, lambda2 = 0.01,
                          beta.naught = rep(0,ncol(data$x)),
                          inner.iter = 100, outer.iter = 100, outer.thresh = 0.0001,
                          gamma = 0.8, step = 1, reset = 10, min.frac= 0.05, verbose = FALSE,
                          groupW = NULL){
  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)

  ## Coming up with other C++ info ##

  num.groups = length(group.length)
  range.group.ind = c(0, cumsum(group.length))
  index = group.length2index(group.length)

  if(is.null(groupW)){
    groupW = sqrt(group.length)
  }

  beta <- beta.naught

  beta.is.zero <- (beta == 0) + 0
  beta.old <- beta

  eta <- rep(0,n)

  intercept <- log(sum(y)) - log(n-sum(y))

  eta = eta + intercept

  junk <- .C("logitNest", X = as.double(as.vector(X)), y = as.integer(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambda1), lambda2 = as.double(lambda2), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), betaZero = as.double(intercept), step = as.double(step), groupW = as.double(groupW))
  if(!is.na(junk$betaZero){
  intercept = junk$betaZero
  beta <- junk$beta
  }
  sol<-list(beta = beta, intercept = intercept )
  return(sol)
}
