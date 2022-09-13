# rCausalMGM

## Requirements

For Windows and Linux users, the rCausalMGM package can be installed directly in R (>=4.0) following the instructions outlined below. MacOS users need to install a `gfortran` compiler in order to build the rCausalMGM package. `gfortran` compilers for a variety of Mac OS versions can be found [here](https://github.com/fxcoudert/gfortran-for-macOS/releases).

## Install and Load rCausalMGM

The rCausalMGM package can be installed directly from this GitHub repository by executing the following code:

```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("tyler-lovelace1/rCausalMGM")
```

To load the rCausalMGM library into your workspace, simply execute

```R
library(rCausalMGM)
```
