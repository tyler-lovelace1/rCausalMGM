#!/bin/bash

R CMD BATCH build.R
# R CMD build .
# R CMD INSTALL rCausalMGM_1.0.tar.gz
R CMD INSTALL rCausalMGM

Rscript rCausalMGM/test/test_MGM.R # > ../test_results/C++debug.txt
