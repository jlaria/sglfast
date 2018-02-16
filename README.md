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

## References

Simon, N., J. Friedman, T. Hastie, and R. Tibshirani (2013). A sparse-group lasso.
*Journal of Computational and Graphical Statistics 22* (2), 231–245.
