# R package `sglfast`

`sglfast` is a fork of **R** package `SGL` (Simon, Friedman, Hastie, Tibshirani 2013), with individual group regularization parameters, and the iterative sparse-group lasso `isgl`, an algorithm to select the optimal regularization parameters of the sparse-group lasso. 

## Installation

The easiest way to install `isgl` in **R** is using `devtools`. 
If `devtools` is not installed, run

```
install.packages('devtools')
```

Then, install `isgl`

```
library(devtools)
install_github("jlaria/sglfast")
library(sglfast)
```

