# rCausalMGM

## Requirements

The Rcpp, RcppArmadillo, and Armadillo packages must be installed. Rcpp and RcppArmadillo can be done in R, but you may need to install Armadillo separately if RcppArmadillo doesn't do it for you.

## Build and Install rCausalMGM

The following steps need to be followed to build and install rCausalMGM as an R library after editing the code.

1. First, we need to generate the `RcppExports.cpp` file that contains the functions that will be exported to R. This only needs to be done when you've defined a new function being exported to R or have changed the output type or parameters of an already defined funtion. In R, execute the following:

```R
Rcpp::compileAtributes("path/to/rCausalMGM")
```

The path can be left blank if you are in the rCausalMGM directory. 

2. Next, you need to build the executable and package it in a way that can be installed in R. To do this, execute the following command in the command line:

```
R CMD build path/to/rCausalMGM
```

3. Finally, it just has to be installed in R so that it recognizes it as a library. To do this, execute the following command in the command line:

```
R CMD INSTALL rCausalMGM_1.0.tar.gz
```

After this has been executed, the R functions in the package can be used by opening R and running:

```R
library(rCausalMGM)
```

For an example, see `test/test_DataSet.R`.