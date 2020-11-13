library(rCausalMGM)

df <- read.table("data/data5.csv", header=T) # Tiny dataset
# df <- read.table("data/data0.txt", header=T) 


# rCausalMGM::runTests(df)
# rCausalMGM::mgm(df)
# g <- rCausalMGM::loadGraph("../cpc_out_censored.txt")
# rCausalMGM::saveGraph(g, "../r_graph_out.txt")
# rCausalMGM::pcStable(df, verbose=TRUE)

g <- rCausalMGM::steps(df, verbose=TRUE)
rCausalMGM::cpcStable(df, initialGraph = g, verbose = TRUE)

# rCausalMGM::pcMax(df, lambda = c(0.05, 0.05, 0.05), verbose = TRUE)

