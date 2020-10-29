library(rCausalMGM)

df <- read.table("data/data5.csv", header=T)
# df <- read.table("data/data0.txt", header=T) # Tiny dataset


rCausalMGM::MGMTest(df)