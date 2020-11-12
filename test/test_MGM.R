library(rCausalMGM)

# df <- read.table("data/data5.csv", header=T) # Tiny dataset
df <- read.table("data/data0.txt", header=T) 


# rCausalMGM::MGMTest(df)
# rCausalMGM::mgm(df)
g <- rCausalMGM::loadGraph("../cpc_out_censored.txt")
rCausalMGM::saveGraph(g, "../r_graph_out.txt")

