library(rCausalMGM)

df <- read.table("data/data0.csv", header=T)

rCausalMGM::MGMTest(df)