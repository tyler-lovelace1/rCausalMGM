library(rCausalMGM)

df <- read.table("data/data5.csv", header=T) # Tiny dataset
# df <- read.table("data/data0.txt", header=T) 


rCausalMGM::MGMTest(df)