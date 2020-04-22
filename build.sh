#!/bin/bash

R CMD BATCH build.R
R CMD build .
R CMD INSTALL rCausalMGM_1.0.tar.gz

Rscript test/test_MGM.R > /mnt/c/Users/maxdu/OneDrive/Benos/Testing/C++debug.txt