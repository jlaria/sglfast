test_error_linear <- function(X, y, model_params){
  diff <- y - (X%*%model_params$beta + model_params$intercept)
  B <- 0.5/length(y)*norm(diff,'F')^2
  return(B)
}

test_error_logit <- function(X, y, model_params){
  model_prediction <-  X%*%model_params$beta + model_params$intercept
  log_sum = sum(log(1+exp(model_prediction)))
  B <- 1/length(y)*(log_sum - t(model_prediction)%*%y)
  return(B)
}
