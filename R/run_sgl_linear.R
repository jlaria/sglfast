run_sgl_linear <- function(data, group.length, thresh = 0.0001, lambda1 = NULL, lambda2=NULL,
                           beta.naught = rep(0,ncol(data$x)),
                           inner.iter = 100, outer.iter = 100,
                           outer.thresh = 0.0001, gamma = 0.8, step = 1, reset = 10,
                           min.frac = 0.05, verbose = FALSE,
                           groupW = NULL){
  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)

  ## Setting up group lasso stuff ##

  num.groups = length(group.length)
  range.group.ind = c(0, cumsum(group.length))
  index = group.length2index(group.length)

  if( is.null(groupW) ){
    groupW = sqrt(group.length)
  }

  ## DONE SETTING UP C STUFF ##

  beta.old <- beta.naught
  beta.is.zero <- (beta.old == 0) + 0
  beta <- rep(0, ncol(X))

  eta <- rep(0,n)

  junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambda1), lambda2 = as.double(lambda2), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset), groupW = as.double(groupW))

  beta <- junk$beta

  sol = list(beta = beta)
  return(sol)

}
