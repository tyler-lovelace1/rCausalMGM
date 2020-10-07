library(rCausalMGM)

df <- read.table("data/data0_smaller.txt", header=T)

rCausalMGM::MGMTest(df)