#!/bin/bash

R CMD BATCH build.R
R CMD build .
R CMD INSTALL rCausalMGM_1.0.tar.gz

Rscript test/test_MGM.R > ../Testing/C++debug.txt