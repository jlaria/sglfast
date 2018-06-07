# R package `sglfast`

`sglfast` is a fork of **R** package `SGL` (Simon et. al. 2013), with individual group regularization parameters, and the iterative sparse-group lasso `isgl`, an algorithm to select the optimal regularization parameters of the sparse-group lasso. 

## Installation

The easiest way to install `isgl` in **R** is using `devtools`. 
If `devtools` is not installed, run


	install.packages('devtools')


Then, install `isgl`

	library(devtools)
	install_github("jlaria/sglfast")
	library(sglfast)

## Usage

Example 1. Estimating `beta` under linear response.

    library(sglfast)
    
    # We create beta="the true coefficient vector" to be used in the simulations.
    beta = 1:5
    
    # We generate the model matrix X with iid columns and rows and the response y
    X = matrix(rnorm(100*400), nrow = 100)
    y = X[,1:5]%*%beta
    
    # We chose the variance of the error such that SNR = 3
    snr = 3
    error = rnorm(100, mean = 0, sd=sqrt(var(y)/snr))
    y = y+error

    # Rows in the training sample
    train.idx = sample(100, 50)
  
    # Group indices for the SGL  
    group_index = rep(1:40, each=10)

    # Input data for the iterative 
    data.train = list(x=X[train.idx,], y=y[train.idx])
    data.validate = list(x=X[-train.idx,], y=y[-train.idx])

    # We run the (unpooled) iterative SGL. For the 2-parameter version use isg_simple()
    isgl.fit = isgl(data.train, data.validate, group_index, type = "linear")
    
    # Best model returned by the iSGL algorithm
    isgl.fit$beta
    isgl.fit$intercept

## References

Simon, N., J. Friedman, T. Hastie, and R. Tibshirani (2013). A sparse-group lasso.
*Journal of Computational and Graphical Statistics 22* (2), 231–245.
