library(rCausalMGM)

# df <- read.table("data/data0_smaller.txt", header=T)
df <- read.table("data/data5.csv", header=T) # Tiny dataset


rCausalMGM::MGMTest(df)